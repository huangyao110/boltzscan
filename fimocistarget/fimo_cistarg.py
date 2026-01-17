"""
使用pyMEMEsuite API构建cisTarget数据库
支持多个.meme文件组成的motif数据库目录
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Union
import glob
from Bio import SeqIO
from pymemesuite.common import MotifFile, Sequence, Background, Alphabet, Array
from pymemesuite.fimo import FIMO
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PyMEMESuiteCisTargetBuilder:
    """
    使用pyMEMEsuite API构建cisTarget数据库
    """
    
    def __init__(self, 
                 fasta_file: str, 
                 motif_dir: str,
                 output_dir: str = './cistarg_db',
                 custom_bg: dict =None,
                 pvalue_thresh: float = 1e-4,
                 max_stored_scores: int = 500000,
                 motif_pseudo: float = 0.1):
        """
        初始化构建器
        
        Parameters:
        -----------
        fasta_file : str
            包含所有区域序列的FASTA文件
        motif_dir : str
            包含.meme文件的目录路径
        output_dir : str
            输出目录
        pvalue_thresh : float
            FIMO的P-value阈值
        max_stored_scores : int
            FIMO存储的最大得分数
        motif_pseudo : float
            Motif伪计数
        """
        self.fasta_file = Path(fasta_file)
        self.motif_dir = Path(motif_dir)
        self.output_dir = Path(output_dir)
        self.pvalue_thresh = pvalue_thresh
        self.max_stored_scores = max_stored_scores
        self.motif_pseudo = motif_pseudo
        if custom_bg is None:
            self.custom_bg = None
        else:
            self.custom_bg = Background(alphabet = Alphabet.dna(), frequencies = Array([custom_bg['A'], custom_bg['C'], 
                                                                                    custom_bg['G'], custom_bg['T']]))
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 检查文件
        if not self.fasta_file.exists():
            raise FileNotFoundError(f"FASTA文件不存在: {self.fasta_file}")
        if not self.motif_dir.exists():
            raise FileNotFoundError(f"Motif目录不存在: {self.motif_dir}")
        
        logger.info(f"初始化完成")
        logger.info(f"FASTA文件: {self.fasta_file}")
        logger.info(f"Motif目录: {self.motif_dir}")
        logger.info(f"P-value阈值: {self.pvalue_thresh}")
    
    def load_motif_files(self) -> List[Path]:
        """
        加载motif目录中所有的.meme文件
        
        Returns:
        --------
        motif_files : List[Path]
            .meme文件路径列表
        """
        # 支持.meme和.txt
        motif_files = list(self.motif_dir.glob('*.meme'))
        motif_files.extend(self.motif_dir.glob('*.txt'))
        
        if len(motif_files) == 0:
            raise ValueError(f"在 {self.motif_dir} 中未找到.meme文件")
        
        motif_files = sorted(motif_files)
        logger.info(f"找到 {len(motif_files)} 个motif文件")
        
        return motif_files
    
    def load_sequences(self) -> List[Sequence]:
        """
        从FASTA文件加载序列
        
        Returns:
        --------
        sequences : List[Sequence]
            pyMEMEsuite的Sequence对象列表
        """
        sequences = [
            Sequence(str(record.seq), name=record.id.encode())
            for record in SeqIO.parse(self.fasta_file, "fasta")
        ]
        logger.info(f"加载了 {len(sequences)} 个序列")
        return sequences
    
    def run_fimo_single_file(self, 
                            motif_file: Path, 
                            sequences: List[Sequence]) -> pd.DataFrame:
        """
        对单个motif文件运行FIMO
        
        Parameters:
        -----------
        motif_file : Path
            单个.meme文件路径
        sequences : List[Sequence]
            序列列表
            
        Returns:
        --------
        results_df : pd.DataFrame
            FIMO结果
        """
        logger.info(f"处理motif文件: {motif_file.name}")
        
        try:
            # 加载motif文件
            with MotifFile(str(motif_file)) as motif_fh:
                motif = motif_fh.read()
            
            fimo = FIMO(both_strands=True)   
            pattern = fimo.score_motif(
                motif,
                sequences,
                self.custom_bg
            )
            
            # 解析结果
            results = []
            for result in pattern.matched_elements:
                results.append({
                    'sequence_name': result.source.accession.decode(),
                    'start': result.start,
                    'stop': result.stop,
                    'strand': result.strand,
                    'score': result.score,
                    'pvalue': result.pvalue,
                    'qvalue': result.qvalue,
                })
            
            if len(results) == 0:
                logger.warning(f"  - 未找到任何匹配")
                return pd.DataFrame()
            
            results_df = pd.DataFrame(results)
            logger.info(f"  - 找到 {len(results_df)} 个匹配")
            
            return results_df
            
        except Exception as e:
            logger.error(f"处理 {motif_file.name} 时出错: {str(e)}")
            return pd.DataFrame()
    
    def run_fimo_all_files(self) -> pd.DataFrame:
        """
        对所有motif文件运行FIMO
        
        Returns:
        --------
        all_results : pd.DataFrame
            合并的FIMO结果
        """
        motif_files = self.load_motif_files()
        sequences = self.load_sequences()
        
        all_results = []
        
        logger.info("=" * 70)
        logger.info("开始运行FIMO扫描...")
        logger.info("=" * 70)
        
        for i, motif_file in enumerate(motif_files, 1):
            logger.info(f"\n[{i}/{len(motif_files)}] 处理文件: {motif_file.name}")
            
            results_df = self.run_fimo_single_file(motif_file, sequences)
            
            if not results_df.empty:
                # 添加来源文件信息
                results_df['motif_id'] = motif_file.stem
                all_results.append(results_df)
        
        if len(all_results) == 0:
            raise ValueError("所有motif文件都没有找到匹配！请检查：\n"
                           "1. P-value阈值是否太严格\n"
                           "2. Motif文件格式是否正确\n"
                           "3. 序列是否有效")
        
        # 合并所有结果
        combined_df = pd.concat(all_results, ignore_index=True)
        
        logger.info("\n" + "=" * 70)
        logger.info("FIMO扫描完成")
        logger.info("=" * 70)
        logger.info(f"总匹配数: {len(combined_df)}")
        logger.info(f"涉及序列数: {combined_df['sequence_name'].nunique()}")
        logger.info(f"涉及motif数: {combined_df['motif_id'].nunique()}")
        
        # 保存原始结果
        output_file = self.output_dir / 'fimo_results_raw.csv.gz'
        combined_df.to_csv(output_file, index=False, compression='gzip')
        logger.info(f"原始结果已保存: {output_file}")
        
        return combined_df
    
    def analyze_multiple_hits(self, fimo_df: pd.DataFrame) -> pd.DataFrame:
        """
        分析多次命中的统计信息
        
        Parameters:
        -----------
        fimo_df : pd.DataFrame
            FIMO结果DataFrame
            
        Returns:
        --------
        hit_stats : pd.DataFrame
            每个区域-motif对的命中统计
        """
        logger.info("\n分析多次命中统计...")
        
        stats = fimo_df.groupby(['sequence_name', 'motif_id']).agg({
            'score': ['count', 'max', 'mean', 'std', 'sum'],
            'pvalue': 'min',
            'start': lambda x: list(x)
        }).reset_index()
        
        stats.columns = ['_'.join(col).strip('_') for col in stats.columns]
        stats = stats.rename(columns={
            'score_count': 'n_hits',
            'score_max': 'max_score',
            'score_mean': 'mean_score',
            'score_std': 'std_score',
            'score_sum': 'sum_score',
            'pvalue_min': 'best_pvalue',
            'start_<lambda>': 'positions'
        })
        
        logger.info(f"平均每个区域-motif对的命中次数: {stats['n_hits'].mean():.2f}")
        logger.info(f"最大命中次数: {stats['n_hits'].max()}")
        
        # 找出命中次数最多的motif
        top_hits = stats.nlargest(10, 'n_hits')[
            ['sequence_name', 'motif_id', 'n_hits', 'max_score', 'mean_score']
        ]
        logger.info("\n命中次数最多的前10个区域-motif对：")
        logger.info(f"\n{top_hits.to_string()}")
        
        # 保存统计信息
        stats_file = self.output_dir / 'hit_statistics.csv'
        stats.to_csv(stats_file, index=False)
        logger.info(f"\n命中统计已保存: {stats_file}")
        
        return stats
    
    def aggregate_scores(self, 
                        fimo_df: pd.DataFrame, 
                        method: str = 'max',
                        top_n: int = 3) -> pd.DataFrame:
        """
        聚合每个区域-motif对的得分
        
        Parameters:
        -----------
        fimo_df : pd.DataFrame
            FIMO结果DataFrame
        method : str
            聚合方法：'max', 'sum', 'mean', 'count', 'top_n_mean', 
                     'weighted_sum', 'max_with_count'
        top_n : int
            用于top_n_mean方法
            
        Returns:
        --------
        score_matrix : pd.DataFrame
            区域×motif的得分矩阵
        """
        logger.info(f"\n使用 '{method}' 方法聚合得分...")
        
        # 添加权重列（基于p-value）
        fimo_df['weight'] = -np.log10(fimo_df['pvalue'] + 1e-300)
        
        if method == 'max':
            score_pivot = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc='max',
                fill_value=0
            )
            
        elif method == 'sum':
            score_pivot = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc='sum',
                fill_value=0
            )
            
        elif method == 'mean':
            score_pivot = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc='mean',
                fill_value=0
            )
            
        elif method == 'count':
            score_pivot = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc='count',
                fill_value=0
            )
            
        elif method == 'top_n_mean':
            def top_n_mean_func(x):
                return x.nlargest(min(top_n, len(x))).mean()
            
            score_pivot = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc=top_n_mean_func,
                fill_value=0
            )
            
        elif method == 'weighted_sum':
            fimo_df['weighted_score'] = fimo_df['score'] * fimo_df['weight']
            score_pivot = fimo_df.pivot_table(
                values='weighted_score',
                index='sequence_name',
                columns='motif_id',
                aggfunc='sum',
                fill_value=0
            )
            
        elif method == 'max_with_count':
            grouped = fimo_df.groupby(['sequence_name', 'motif_id'])
            
            def max_with_count(group):
                max_score = group['score'].max()
                count = len(group)
                return max_score * np.log(count + 1)
            
            agg_scores = grouped.apply(max_with_count)
            score_pivot = agg_scores.unstack(fill_value=0)
            
        else:
            raise ValueError(f"不支持的聚合方法: {method}")
        
        logger.info(f"得分矩阵维度: {score_pivot.shape}")
        logger.info(f"  - {len(score_pivot)} 个区域")
        logger.info(f"  - {len(score_pivot.columns)} 个motifs")
        
        return score_pivot
    
    def create_ranking_db(self, score_matrix: pd.DataFrame) -> pd.DataFrame:
        """
        从得分矩阵创建排序数据库
        
        Parameters:
        -----------
        score_matrix : pd.DataFrame
            区域×motif的得分矩阵
            
        Returns:
        --------
        ranking_db : pd.DataFrame
            区域×motif的排序矩阵（1=最高得分）
        """
        logger.info("\n创建排序数据库...")
        
        # 对每个motif（列）进行排序
        ranking_db = score_matrix.rank(ascending=False, method='min')
        
        logger.info(f"排序数据库维度: {ranking_db.shape}")
        
        return ranking_db
    
    def save_databases(self, 
                      score_matrix: pd.DataFrame, 
                      ranking_db: pd.DataFrame):
        """
        保存得分矩阵和排序数据库
        
        Parameters:
        -----------
        score_matrix : pd.DataFrame
            得分矩阵
        ranking_db : pd.DataFrame
            排序数据库
        """
        logger.info("\n保存数据库文件...")
        
        # 保存为feather格式（快速）
        score_feather = self.output_dir / 'score_matrix.feather'
        ranking_feather = self.output_dir / 'ranking_db.feather'
        
        score_matrix.reset_index().to_feather(score_feather)
        ranking_db.reset_index().to_feather(ranking_feather)
        
        logger.info(f"  - 得分矩阵 (feather): {score_feather}")
        logger.info(f"  - 排序数据库 (feather): {ranking_feather}")
        
        # 同时保存CSV格式（可读性）
        score_csv = self.output_dir / 'score_matrix.csv.gz'
        ranking_csv = self.output_dir / 'ranking_db.csv.gz'
        
        score_matrix.to_csv(score_csv, compression='gzip')
        ranking_db.to_csv(ranking_csv, compression='gzip')
        
        logger.info(f"  - 得分矩阵 (csv.gz): {score_csv}")
        logger.info(f"  - 排序数据库 (csv.gz): {ranking_csv}")
    
    def build_database(self, 
                      aggregation_method: str = 'max',
                      analyze_hits: bool = True,
                      top_n: int = 3) -> tuple:
        """
        完整的数据库构建流程
        
        Parameters:
        -----------
        aggregation_method : str
            得分聚合方法
        analyze_hits : bool
            是否分析多次命中统计
        top_n : int
            用于top_n_mean方法
            
        Returns:
        --------
        score_matrix, ranking_db : tuple
            得分矩阵和排序数据库
        """
        logger.info("\n" + "=" * 70)
        logger.info("开始构建cisTarget数据库")
        logger.info("=" * 70)
        
        # 步骤1: 运行FIMO
        fimo_df = self.run_fimo_all_files()
        
        # 步骤2: 分析多次命中
        if analyze_hits:
            self.analyze_multiple_hits(fimo_df)
        
        # 步骤3: 聚合得分
        score_matrix = self.aggregate_scores(
            fimo_df, 
            method=aggregation_method,
            top_n=top_n
        )
        
        # 步骤4: 创建排序数据库
        ranking_db = self.create_ranking_db(score_matrix)
        
        # 步骤5: 保存数据库
        self.save_databases(score_matrix, ranking_db)
        
        logger.info("\n" + "=" * 70)
        logger.info("数据库构建完成！")
        logger.info("=" * 70)
        
        return score_matrix, ranking_db



# 命令行接口
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='使用pyMEMEsuite构建cisTarget数据库'
    )
    parser.add_argument('--fasta', required=True, 
                       help='区域序列FASTA文件')
    parser.add_argument('--motif-dir', required=True, 
                       help='包含.meme文件的目录')
    parser.add_argument('--output', default='./cistarg_db', 
                       help='输出目录')
    parser.add_argument('--pvalue', type=float, default=1e-4, 
                       help='FIMO P-value阈值')
    parser.add_argument('--method', default='max', 
                       choices=['max', 'sum', 'mean', 'top_n_mean', 
                               'max_with_count', 'weighted_sum'],
                       help='得分聚合方法')
    parser.add_argument('--top-n', type=int, default=3,
                       help='top_n_mean方法的N值')
    parser.add_argument('--no-analyze', action='store_true',
                       help='不进行多次命中分析（加快速度）')
    
    args = parser.parse_args()
    
    # 构建数据库
    builder = PyMEMESuiteCisTargetBuilder(
        fasta_file=args.fasta,
        motif_dir=args.motif_dir,
        output_dir=args.output,
        pvalue_thresh=args.pvalue
    )
    
    score_matrix, ranking_db = builder.build_database(
        aggregation_method=args.method,
        analyze_hits=not args.no_analyze,
        top_n=args.top_n
    )
    
    print("\n完成！数据库文件已保存到:", args.output)
