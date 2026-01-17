import pandas as pd
import numpy as np
from scipy import stats
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import logging
from Bio import SeqIO
# 如果需要多重检验校正，需要安装statsmodels
# pip install statsmodels
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GenomeBackgroundComparator:
    """
    基因组背景比较分析器
    """
    
    # ... (原有方法 __init__, _validate_data, aggregate_scores_per_sequence, etc. 保持不变) ...
    def __init__(self, 
                 background_fimo_results: pd.DataFrame,
                 target_fimo_results: pd.DataFrame,
                 score_aggregation: str = 'max',
                 gc_correction: bool = True):
        """
        初始化分析器
        
        Parameters:
        -----------
        background_fimo_results : pd.DataFrame
            背景启动子的FIMO扫描结果
            必须包含列：'sequence_name', 'motif_id', 'score', 'pvalue'
        target_fimo_results : pd.DataFrame
            目标启动子的FIMO扫描结果
            score_aggregation : str
            得分聚合方法：'max', 'sum', 'mean', 'top_n_mean'
        gc_correction : bool
            是否进行GC含量校正 (注意：在motif-centric分析中，GC校正较为复杂，通常使用全局背景)
        """
        self.background_fimo = background_fimo_results
        self.target_fimo = target_fimo_results
        self.score_aggregation = score_aggregation
        self.gc_correction = gc_correction
        
        # 验证数据
        self._validate_data()
        
        logger.info("基因组背景比较分析器初始化完成")
        logger.info(f"  背景启动子数: {self.background_fimo['sequence_name'].nunique()}")
        logger.info(f"  目标启动子数: {self.target_fimo['sequence_name'].nunique()}")
        logger.info(f"  Motif数: {self.background_fimo['motif_id'].nunique()}")

        # 预计算聚合得分，避免重复计算
        self.background_scores_agg = self.aggregate_scores_per_sequence(self.background_fimo)
        self.target_scores_agg = self.aggregate_scores_per_sequence(self.target_fimo)

    
    def _validate_data(self):
        """验证输入数据格式"""
        required_cols = ['sequence_name', 'motif_id', 'score', 'pvalue']
        
        for col in required_cols:
            if col not in self.background_fimo.columns:
                raise ValueError(f"背景FIMO结果缺少列: {col}")
            if col not in self.target_fimo.columns:
                raise ValueError(f"目标FIMO结果缺少列: {col}")
    
    def aggregate_scores_per_sequence(self, 
                                     fimo_df: pd.DataFrame,
                                     top_n: int = 3) -> pd.DataFrame:
        """
        聚合每个序列-motif对的得分（处理多次命中）
        
        Parameters:
        -----------
        fimo_df : pd.DataFrame
            FIMO结果
        top_n : int
            用于top_n_mean方法
            
        Returns:
        --------
        aggregated : pd.DataFrame
            聚合后的得分矩阵（sequence × motif）
        """
        logger.info(f"使用 {self.score_aggregation} 方法聚合得分...")
        
        if self.score_aggregation == 'max':
            agg_df = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc='max',
                fill_value=0
            )
        
        elif self.score_aggregation == 'sum':
            agg_df = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc='sum',
                fill_value=0
            )
        
        elif self.score_aggregation == 'mean':
            agg_df = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc='mean',
                fill_value=0
            )
        
        elif self.score_aggregation == 'top_n_mean':
            def top_n_mean(x):
                return x.nlargest(min(top_n, len(x))).mean()
            
            agg_df = fimo_df.pivot_table(
                values='score',
                index='sequence_name',
                columns='motif_id',
                aggfunc=top_n_mean,
                fill_value=0
            )
        
        else:
            raise ValueError(f"不支持的聚合方法: {self.score_aggregation}")
        
        return agg_df

    # ... (其他辅助方法如 calculate_gc_content, match_gc_background 等保持不变) ...
    def calculate_gc_content(self, sequences: Dict[str, str]) -> pd.Series:
        gc_dict = {}
        for name, seq in sequences.items():
            seq_upper = seq.upper()
            gc_count = seq_upper.count('G') + seq_upper.count('C')
            gc_dict[name] = gc_count / len(seq) if len(seq) > 0 else 0
        return pd.Series(gc_dict)

    def match_gc_background(self, target_gc: float, background_gc: pd.Series, tolerance: float = 0.1) -> List[str]:
        lower, upper = target_gc - tolerance, target_gc + tolerance
        matched = background_gc[(background_gc >= lower) & (background_gc <= upper)]
        logger.info(f"GC匹配筛选: 目标GC={target_gc:.3f}, 范围=[{lower:.3f}, {upper:.3f}], 保留 {len(matched)}/{len(background_gc)} 背景序列")
        return matched.index.tolist()

    def calculate_background_statistics(self, background_scores: pd.DataFrame) -> pd.DataFrame:
        logger.info("计算背景统计分布...")
        stats_list = []
        for motif in background_scores.columns:
            scores = background_scores[motif].values
            scores_nonzero = scores[scores > 0]
            stats_dict = {
                'motif': motif, 'n_total': len(scores), 'n_detected': len(scores_nonzero),
                'detection_rate': len(scores_nonzero) / len(scores) if len(scores) > 0 else 0,
                'mean': np.mean(scores), 'mean_nonzero': np.mean(scores_nonzero) if len(scores_nonzero) > 0 else 0,
                'median': np.median(scores), 'median_nonzero': np.median(scores_nonzero) if len(scores_nonzero) > 0 else 0,
                'std': np.std(scores), 'std_nonzero': np.std(scores_nonzero) if len(scores_nonzero) > 0 else 0,
                'min': np.min(scores), 'max': np.max(scores), 'q25': np.percentile(scores, 25),
                'q75': np.percentile(scores, 75), 'q90': np.percentile(scores, 90),
                'q95': np.percentile(scores, 95), 'q99': np.percentile(scores, 99)
            }
            stats_list.append(stats_dict)
        return pd.DataFrame(stats_list)

    def calculate_zscore(self, target_score: float, background_mean: float, background_std: float) -> float:
        return (target_score - background_mean) / background_std if background_std > 0 else 0.0

    def calculate_percentile_rank(self, target_score: float, background_scores: np.ndarray) -> float:
        return stats.percentileofscore(background_scores, target_score) if len(background_scores) > 0 else 50.0

    def calculate_empirical_pvalue(self, target_score: float, background_scores: np.ndarray, alternative: str = 'greater') -> float:
        n = len(background_scores)
        if n == 0: return 1.0
        if alternative == 'greater':
            count = np.sum(background_scores >= target_score)
        else:
            count = np.sum(background_scores <= target_score)
        return (count + 1) / (n + 1)

    def compare_single_promoter(self,
                               target_name: str,
                               background_sequences: Dict[str, str] = None,
                               gc_tolerance: float = 0.1,
                               use_nonzero_only: bool = False) -> pd.DataFrame:
        """
        比较单个目标启动子与基因组背景 (Promoter-centric analysis)
        """
        logger.info(f"\n分析目标启动子: {target_name}")
        logger.info("=" * 70)
        
        background_scores = self.background_scores_agg.copy()
        
        if self.gc_correction and background_sequences is not None:
            logger.info("\n应用GC含量匹配...")
            background_gc = self.calculate_gc_content(background_sequences)
            if target_name in background_sequences:
                target_seq = background_sequences[target_name]
                target_gc = (target_seq.upper().count('G') + target_seq.upper().count('C')) / len(target_seq)
                matched_names = self.match_gc_background(target_gc, background_gc, gc_tolerance)
                matched_names = [name for name in matched_names if name != target_name]
                if len(matched_names) > 0:
                    background_scores = background_scores.loc[background_scores.index.intersection(matched_names)]
                else:
                    logger.warning("GC匹配后没有剩余的背景序列，将使用全部背景。")
            else:
                logger.warning(f"目标序列 {target_name} 不在序列字典中, 跳过GC匹配")
        
        bg_stats = self.calculate_background_statistics(background_scores)
        
        if target_name not in self.target_scores_agg.index:
            raise ValueError(f"目标启动子 {target_name} 不在FIMO结果中")
        
        target_row = self.target_scores_agg.loc[target_name]
        
        results = []
        for motif in bg_stats['motif']:
            target_score = target_row.get(motif, 0)
            bg_info = bg_stats[bg_stats['motif'] == motif].iloc[0]
            
            bg_mean = bg_info['mean_nonzero'] if use_nonzero_only else bg_info['mean']
            bg_median = bg_info['median_nonzero'] if use_nonzero_only else bg_info['median']
            bg_std = bg_info['std_nonzero'] if use_nonzero_only else bg_info['std']
            
            zscore = self.calculate_zscore(target_score, bg_mean, bg_std)
            bg_scores_array = background_scores[motif].values
            percentile = self.calculate_percentile_rank(target_score, bg_scores_array)
            emp_pvalue = self.calculate_empirical_pvalue(target_score, bg_scores_array, 'greater')
            
            results.append({
                'motif': motif, 'target_score': target_score, 'bg_mean': bg_mean, 'bg_std': bg_std,
                'bg_median': bg_info['median'], 'bg_q95': bg_info['q95'], 'bg_q99': bg_info['q99'],
                'zscore': zscore, 'percentile': percentile, 'empirical_pvalue': emp_pvalue,
                'bg_detection_rate': bg_info['detection_rate'],
                'significant': (zscore > 2.0) and (target_score > bg_mean),
                'highly_significant': (zscore > 3.0) and (percentile > 99)
            })
        
        results_df = pd.DataFrame(results).sort_values('zscore', ascending=False)
        logger.info(f"\n分析完成: 显著富集Motifs (Z>2): {results_df['significant'].sum()}, "
                    f"高度显著 (Z>3, perc>99): {results_df['highly_significant'].sum()}")
        return results_df

    # ... (原有 get_significant_motifs, visualize_comparison 方法保持不变) ...
    def get_significant_motifs(self,
                            comparison_results: pd.DataFrame,
                            zscore_threshold: float = 2.0,
                            min_target_score: float = 8.0) -> pd.DataFrame:
        """
        筛选显著的motifs
        
        Parameters:
        -----------
        comparison_results : pd.DataFrame
            比较结果
        zscore_threshold : float
            Z-score阈值
        min_target_score : float
            最小目标得分
            
        Returns:
        --------
        significant : pd.DataFrame
            显著的motifs
        """
        significant = comparison_results[
            (comparison_results['zscore'] > zscore_threshold) &
            (comparison_results['target_score'] > min_target_score)
        ].copy()
        
        logger.info(f"\nSignificant motifs filtering:")
        logger.info(f"  Z-score > {zscore_threshold}")
        logger.info(f"  Target score > {min_target_score}")
        logger.info(f"  Result: {len(significant)} motifs")
        
        return significant

    def visualize_comparison(self,
                            comparison_results: pd.DataFrame,
                            top_n: int = 20,
                            output_file: str = 'motif_comparison.png'):
        """
        可视化比较结果
        
        Parameters:
        -----------
        comparison_results : pd.DataFrame
            比较结果
        top_n : int
            显示前N个motifs
        output_file : str
            输出文件名
        """
        # 取前N个Z-score最高的motifs
        top_motifs = comparison_results.nlargest(top_n, 'zscore')
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 图1：Z-score条形图
        ax1 = axes[0, 0]
        colors = ['red' if z > 3 else 'orange' if z > 2 else 'gray' 
                for z in top_motifs['zscore']]
        ax1.barh(range(len(top_motifs)), top_motifs['zscore'], color=colors)
        ax1.set_yticks(range(len(top_motifs)))
        ax1.set_yticklabels(top_motifs['motif'], fontsize=10)
        ax1.axvline(x=2, color='orange', linestyle='--', label='Z=2 (P<0.05)')
        ax1.axvline(x=3, color='red', linestyle='--', label='Z=3 (P<0.001)')
        ax1.set_xlabel('Z-score', fontsize=12)
        ax1.set_title('Motif Enrichment Z-scores', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(axis='x', alpha=0.3)
        ax1.invert_yaxis()
        
        # 图2：目标 vs 背景得分散点图
        ax2 = axes[0, 1]
        scatter = ax2.scatter(top_motifs['bg_mean'], top_motifs['target_score'],
                            c=top_motifs['zscore'], cmap='RdYlGn', s=100, alpha=0.7,
                            edgecolors='black', linewidth=0.5)
        
        # 添加对角线
        max_val = max(top_motifs['bg_mean'].max(), top_motifs['target_score'].max())
        ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, label='y=x')
        
        ax2.set_xlabel('Background Mean Score', fontsize=12)
        ax2.set_ylabel('Target Score', fontsize=12)
        ax2.set_title('Target vs Background Scores', fontsize=14, fontweight='bold')
        plt.colorbar(scatter, ax=ax2, label='Z-score')
        ax2.legend()
        ax2.grid(alpha=0.3)
        
        # 图3：百分位排名
        ax3 = axes[1, 0]
        colors = ['red' if p > 99 else 'orange' if p > 95 else 'gray' 
                for p in top_motifs['percentile']]
        ax3.barh(range(len(top_motifs)), top_motifs['percentile'], color=colors)
        ax3.set_yticks(range(len(top_motifs)))
        ax3.set_yticklabels(top_motifs['motif'], fontsize=10)
        ax3.axvline(x=95, color='orange', linestyle='--', label='95th percentile')
        ax3.axvline(x=99, color='red', linestyle='--', label='99th percentile')
        ax3.set_xlabel('Percentile Rank', fontsize=12)
        ax3.set_xlim(0, 100)
        ax3.set_title('Percentile Rankings', fontsize=14, fontweight='bold')
        ax3.legend()
        ax3.grid(axis='x', alpha=0.3)
        ax3.invert_yaxis()
        
        # 图4：经验P值（-log10）
        ax4 = axes[1, 1]
        neg_log_p = -np.log10(top_motifs['empirical_pvalue'] + 1e-10)
        colors = ['red' if p < 0.001 else 'orange' if p < 0.05 else 'gray' 
                for p in top_motifs['empirical_pvalue']]
        ax4.barh(range(len(top_motifs)), neg_log_p, color=colors)
        ax4.set_yticks(range(len(top_motifs)))
        ax4.set_yticklabels(top_motifs['motif'], fontsize=10)
        ax4.axvline(x=-np.log10(0.05), color='orange', linestyle='--', label='P=0.05')
        ax4.axvline(x=-np.log10(0.001), color='red', linestyle='--', label='P=0.001')
        ax4.set_xlabel('-log10(Empirical P-value)', fontsize=12)
        ax4.set_title('Statistical Significance', fontsize=14, fontweight='bold')
        ax4.legend()
        ax4.grid(axis='x', alpha=0.3)
        ax4.invert_yaxis()
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"\n Image had been saved: {output_file}")
        plt.close()

    # ==============================================================================
    # 新增方法：Motif为中心的分析
    # ==============================================================================
    
    def analyze_motif_enrichment_across_promoters(self,
                                                 motif_id: str,
                                                 correct_pvalues: bool = True) -> pd.DataFrame:
        """
        分析单个Motif在所有目标启动子中的富集情况 (Motif-centric analysis)
        
        Parameters:
        -----------
        motif_id : str
            要分析的目标Motif的ID
        correct_pvalues : bool
            是否对经验P值进行多重检验校正 (Benjamini-Hochberg)
            
        Returns:
        --------
        results_df : pd.DataFrame
            包含每个目标启动子相对于背景的统计信息的DataFrame, 按z-score排序
        """
        logger.info(f"\n分析Motif: {motif_id}")
        logger.info("=" * 70)

        # 检查Motif是否存在
        if motif_id not in self.background_scores_agg.columns:
            raise ValueError(f"Motif '{motif_id}' 不在背景FIMO结果中")
        if motif_id not in self.target_scores_agg.columns:
            logger.warning(f"Motif '{motif_id}' 不在目标FIMO结果中, 将视作所有目标得分为0")
            
        # 提取该Motif在背景中的得分分布
        background_motif_scores = self.background_scores_agg[motif_id].values
        bg_mean = np.mean(background_motif_scores)
        bg_std = np.std(background_motif_scores)
        logger.info(f"背景分布 (n={len(background_motif_scores)}): Mean={bg_mean:.2f}, Std={bg_std:.2f}")

        # 获取该Motif在所有目标启动子上的得分
        if motif_id in self.target_scores_agg.columns:
            target_motif_scores = self.target_scores_agg[motif_id]
        else:
            # 如果目标集中完全没有这个motif，创建一个全为0的Series
            target_motif_scores = pd.Series(0, index=self.target_scores_agg.index, name=motif_id)
        
        results = []
        for promoter_name, target_score in target_motif_scores.items():
            zscore = self.calculate_zscore(target_score, bg_mean, bg_std)
            percentile = self.calculate_percentile_rank(target_score, background_motif_scores)
            emp_pvalue = self.calculate_empirical_pvalue(target_score, background_motif_scores, 'greater')
            
            results.append({
                'promoter': promoter_name,
                'target_score': target_score,
                'zscore': zscore,
                'percentile': percentile,
                'empirical_pvalue': emp_pvalue
            })
            
        results_df = pd.DataFrame(results)
        
        # 多重检验校正
        if correct_pvalues and not results_df.empty:
            pvals = results_df['empirical_pvalue']
            reject, pvals_corrected, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
            results_df['p_adjusted'] = pvals_corrected
            results_df['is_significant_adj'] = reject

        # 按z-score排序
        results_df = results_df.sort_values('zscore', ascending=False).reset_index(drop=True)
        
        logger.info(f"\n分析完成: {len(results_df)} 个目标启动子已排序")
        if correct_pvalues and not results_df.empty:
            logger.info(f"  在FDR<0.05水平下显著的启动子数: {results_df['is_significant_adj'].sum()}")
            
        return results_df

    def visualize_motif_enrichment_for_motif(self,
                                             motif_enrichment_results: pd.DataFrame,
                                             motif_id: str,
                                             top_n: int = 20,
                                             output_file: str = None):
        """
        可视化单个Motif在Top N个启动子上的富集情况
        """
        top_promoters = motif_enrichment_results.head(top_n)
        
        plt.figure(figsize=(10, top_n * 0.4))
        
        colors = ['red' if z > 3 else 'orange' if z > 2 else 'skyblue' for z in top_promoters['zscore']]
        
        plt.barh(top_promoters['promoter'], top_promoters['zscore'], color=colors)
        plt.xlabel('Z-score (Enrichment over background)')
        plt.ylabel('Promoter Name')
        plt.title(f'Top {top_n} Promoters Enriched for Motif: {motif_id}', fontsize=16, fontweight='bold')
        
        plt.axvline(2, color='orange', linestyle='--', lw=1, label='Z=2 (Significant)')
        plt.axvline(3, color='red', linestyle='--', lw=1, label='Z=3 (Highly Significant)')
        plt.legend()
        plt.grid(axis='x', linestyle='--', alpha=0.6)
        plt.gca().invert_yaxis() # 将Z-score最高的放在顶部
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Motif富集图像已保存: {output_file}")
            plt.close()
        else:
            plt.show()



import pandas as pd

def remove_redundant_intervals_per_motif(df, overlap_threshold=0.0):
    """
    针对 FIMO 结果去重 (支持重叠比例阈值)：
    逻辑：在 同一个序列(sequence_name) 且 同一个Motif(motif_id) 内部，
          保留分数最高的；如果较低分数的区间与高分区间重叠比例超过 overlap_threshold，则剔除。
    
    参数:
    df: FIMO结果 DataFrame
    overlap_threshold: 0.0 ~ 1.0 (例如 0.5 代表重叠超过50%才剔除，0.0代表只要碰上就剔除)
    """
    # 1. 备份数据
    data = df.copy()

    # 2. 修正坐标：确保 Start < End，并计算长度
    # FIMO 坐标通常是包含的 (inclusive)，所以长度建议 +1 (视具体版本而定，+1更稳妥)
    data['real_start'] = data[['start', 'stop']].min(axis=1)
    data['real_stop'] = data[['start', 'stop']].max(axis=1)
    data['length'] = data['real_stop'] - data['real_start'] + 1 # 计算自身长度

    # 3. 全局按分数降序排序 (优先级核心)
    data = data.sort_values(by=['score'], ascending=False)

    keep_indices = []

    # 4. 按序列和Motif分组处理
    grouped = data.groupby(['sequence_name', 'motif_id'])

    for name, group in grouped:
        accepted_intervals = [] # 存储已保留的区间信息：(start, stop, length)
        
        for idx, row in group.iterrows():
            curr_start = row['real_start']
            curr_stop = row['real_stop']
            curr_len = row['length']
            
            is_redundant = False
            
            # 检查与“本组已保留区间”的重叠情况
            for (acc_start, acc_stop) in accepted_intervals:
                # 1. 计算重叠部分的物理坐标
                # 重叠起点 = 两个起点中的较大值
                # 重叠终点 = 两个终点中的较小值
                overlap_start = max(curr_start, acc_start)
                overlap_end = min(curr_stop, acc_stop)
                
                # 2. 判断是否存在重叠
                if overlap_start < overlap_end:
                    # 计算重叠长度 (+1 是因为坐标是inclusive的)
                    overlap_len = overlap_end - overlap_start + 1
                    
                    # 3. 计算重叠比例
                    # 这里计算的是：重叠部分占“当前这个较低分区间”的比例
                    # 也可以改成占 min(curr_len, acc_len)
                    ratio = overlap_len / curr_len
                    
                    # 4. 如果比例超过阈值，标记为冗余
                    # 注意：如果 threshold 是 0，只要 overlap_len > 0 就会剔除
                    if ratio > overlap_threshold:
                        is_redundant = True
                        break
            
            if not is_redundant:
                accepted_intervals.append((curr_start, curr_stop))
                keep_indices.append(idx)

    # 5. 返回结果，移除临时辅助列
    result = df.loc[keep_indices].sort_values(by=['sequence_name', 'start'])
    return result