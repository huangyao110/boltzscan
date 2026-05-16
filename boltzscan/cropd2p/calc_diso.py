#!/usr/bin/env python3
"""
Protein Disorder Prediction Tool using IDP-EDL Meta Predictor

This module provides functionality to predict protein disorder regions using
a meta predictor that combines IDP, SDR, and LDR task-specific predictors.
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import numpy as np

import torch
import torch.nn as nn
from transformers import T5Tokenizer

from .idpedl.generator import (
    SDRTaskSpecificPredictor,
    LDRTaskSpecificPredictor,
    IDPTaskSpecificPredictor,
    MetaPredictor
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)



class ProteinDisorderPredictor:
    """Protein disorder prediction using IDP-EDL meta predictor."""
    
    def __init__(self, model_dir: str = "./cropd2p/cropdock/idpedl/model", device: Optional[str] = None):
        """
        Initialize the protein disorder predictor.
        
        Args:
            model_dir: Directory containing model weights
            device: Device to run predictions on ('cuda', 'cpu', or None for auto-detection)
        """
        self.model_dir = Path(model_dir)
        self.device = self._get_device(device)
        self.tokenizer = None
        self.meta_predictor = None
        
        logger.info(f"Using device: {self.device}")
        
    def _get_device(self, device: Optional[str] = None) -> torch.device:
        """Get the appropriate device for computation."""
        if device is None:
            return torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        return torch.device(device)
    
    def _load_model_weights(self, model: nn.Module, weights_path: Path) -> nn.Module:
        """
        Load model weights from file.
        
        Args:
            model: PyTorch model to load weights into
            weights_path: Path to the weights file
            
        Returns:
            Model with loaded weights
            
        Raises:
            FileNotFoundError: If weights file doesn't exist
            RuntimeError: If there's an error loading weights
        """
        if not weights_path.exists():
            raise FileNotFoundError(f"Model weights not found: {weights_path}")
        
        try:
            logger.info(f"Loading weights from: {weights_path}")
            params = torch.load(weights_path, weights_only=True, map_location=self.device)
            
            for param_name, param in model.named_parameters():
                if param_name in params:
                    param.data = params[param_name].data
                    
            logger.info(f"Successfully loaded weights for {len(params)} parameters")
            return model
            
        except Exception as e:
            raise RuntimeError(f"Error loading model weights from {weights_path}: {str(e)}")
    
    def load_predictors(self) -> None:
        """Load all component predictors and meta predictor."""
        logger.info("Loading component predictors...")
        
        try:
            # Load component predictors
            generic_predictor, self.tokenizer = IDPTaskSpecificPredictor()
            sdr_predictor, _ = SDRTaskSpecificPredictor()
            ldr_predictor, _ = LDRTaskSpecificPredictor()
            
            # Load pre-trained weights
            generic_predictor = self._load_model_weights(
                generic_predictor, self.model_dir / "generic_predictor.pth"
            )
            sdr_predictor = self._load_model_weights(
                sdr_predictor, self.model_dir / "sdr_predictor.pth"
            )
            ldr_predictor = self._load_model_weights(
                ldr_predictor, self.model_dir / "ldr_predictor.pth"
            )
            
            # Create and configure meta predictor
            self.meta_predictor = MetaPredictor(generic_predictor, sdr_predictor, ldr_predictor)
            self.meta_predictor = self.meta_predictor.to(self.device)
            
            # Load meta predictor weights
            meta_weights_path = self.model_dir / "meta_predictor.pth"
            self._load_model_weights(self.meta_predictor, meta_weights_path)
            
            logger.info("All predictors loaded successfully")
            
        except Exception as e:
            logger.error(f"Error loading predictors: {str(e)}")
            raise RuntimeError(f"Failed to load predictors: {str(e)}")
    
    def preprocess_sequence(self, sequence: str) -> str:
        """
        Preprocess protein sequence for prediction.
        
        Args:
            sequence: Raw protein sequence string
            
        Returns:
            Preprocessed sequence string
        """
        # Remove any whitespace and convert to uppercase
        sequence = " ".join(sequence)
        return sequence
    
    def predict_disorder(self, sequence: str, save_results: bool = False, save_path: Optional[str] = None) -> Dict[str, Union[List[int], torch.Tensor]]:
        """
        Predict protein disorder regions using the meta predictor.
        
        Args:
            sequence: Raw protein sequence string
            save_results: Whether to save results to file
            save_path: Path to save results (required if save_results is True)
            
        Returns:
            Dictionary containing predictions and probabilities
            
        Raises:
            RuntimeError: If predictors are not loaded
            ValueError: If save_results is True but save_path is None
        """
        if self.meta_predictor is None or self.tokenizer is None:
            raise RuntimeError("Predictors not loaded. Call load_predictors() first.")
        
        if save_results and save_path is None:
            raise ValueError("save_path must be provided when save_results is True")
        
        # Preprocess the sequence
        processed_seq = self.preprocess_sequence(sequence)
        
        # Tokenize the sequence
        inputs = self.tokenizer(processed_seq, return_tensors="pt", max_length=1024, padding=False, truncation=True).to(self.device)

        with torch.no_grad():
            outputs = self.meta_predictor(input_ids=inputs["input_ids"], attention_mask=inputs["attention_mask"], labels=None)

        # Extract results
        logits = outputs.logits
        prediction = logits.argmax(dim=-1).tolist()[0][1:]
        assert len(prediction) == len(sequence), f'{len(prediction)}!={len(sequence)}'    
        probabilities = torch.softmax(logits, dim=-1)
        
        results = {
            "predictions": prediction,
            "probabilities": probabilities,
            "sequence": processed_seq
        }
        
        if save_results:
            self._save_results(results, save_path)
            
        return results
    
    def pred_from_fasta(self, fasta_path: str, save_dir: Optional[str] = None) -> Dict[str, Union[List[int], torch.Tensor]]:
        """
        Predict protein disorder regions from a FASTA file.

        Args:
            fasta_path: Path to FASTA file
            save_results: Whether to save results to file
            save_path: Path to save results (required if save_results is True)

        """
        fasta_path = Path(fasta_path)
        save_dir = Path(save_dir) if save_dir else None
        from Bio import SeqIO
        seq_dct = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
        for seq_id, seq in seq_dct.items():
            seq = str(seq.seq)
            save_path = save_dir / f"{seq_id}.npz" if save_dir else None
            _ = self.predict_disorder(sequence=seq, save_results=True, save_path=save_path)
            
        logger.info(f"Predictions saved to: {save_dir}")
    
    def _save_results(self, data: Dict, filename: str) -> None:
        """
        Save prediction results to NPZ file.
        
        Args:
            data: Dictionary containing predictions and probabilities
            filename: Output filename
        """
        # Convert tensors to numpy arrays for saving
        save_data = {}
        for key, value in data.items():
            if isinstance(value, torch.Tensor):
                save_data[key] = value.cpu().numpy()
            else:
                save_data[key] = value
        
        np.savez(filename, **save_data)
        logger.info(f"Results saved to: {filename}")


def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description="Protein Disorder Prediction Tool using IDP-EDL Meta Predictor"
    )
    parser.add_argument(
        "-f",
        "--fasta_file",
        help="Protein sequence string or path to FASTA file"
    )
    parser.add_argument(
        "-md",
        "--model-dir",
        default="./IDP-EDL/model",
        help="Directory containing model weights (default: ./IDP-EDL/model)"
    )
    parser.add_argument(
        "-d",
        "--device",
        choices=["cuda", "cpu"],
        help="Device to run predictions on (default: auto-detect)"
    )
    parser.add_argument(
        "-o"
        "--output",
        help="Output file path for saving results (required if --save-results is used)"
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        # Initialize predictor
        predictor = ProteinDisorderPredictor(model_dir=args.model_dir, device=args.device)
        
        # Load predictors
        logger.info("Loading predictors...")
        predictor.load_predictors()
        
        # Get sequence
        f = Path(args.fasta_file)
        
        # Make prediction
        logger.info("Making prediction...")
        results = predictor.pred_from_fasta(
            fasta_path=f,
            save_path=args.output
        )
        
        # Print results summary
        logger.info("Prediction completed successfully!")
        logger.info(f"Sequence length: {len(results['sequence'])}")
        logger.info(f"Number of predictions: {len(results['predictions'])}")
        
        if not args.save_results:
            logger.info("Use --save-results and --output to save detailed results")
            
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()