import pandas as pd

class Filter:

    def __init__(self, base_dataframe, pval_range, min_imputation_info_score, remove_strand_ambig_snps) -> None:
        self.base_dataframe = base_dataframe
        self.pval_range = pval_range
        self.min_imputation_info_score = min_imputation_info_score
        self.remove_strand_ambig_snps = remove_strand_ambig_snps
        self.filter()

    def filter(self):
        summary_stats = self.base_dataframe
        snp_counter = summary_stats.shape[0]
        print(f'{snp_counter:,} SNPs in complete dataframe')

        if self.pval_range[0] == None:
            summary_stats = summary_stats[summary_stats.P < self.pval_range[-1]].reset_index(drop=True)
            print(f'-{snp_counter - summary_stats.shape[0]:,} SNPs are out of the {self.pval_range} p-value bounds ({summary_stats.shape[0]:,} SNPs remain)')
            snp_counter = summary_stats.shape[0]
        else:        
            min_p_filt = summary_stats.P > self.pval_range[0] # P-value filter
            max_p_filt = summary_stats.P < self.pval_range[-1]
            summary_stats = summary_stats[min_p_filt & max_p_filt].reset_index(drop=True)
            print(f'-{snp_counter - summary_stats.shape[0]:,} SNPs are out of the {self.pval_range} p-value bounds ({summary_stats.shape[0]:,} SNPs remain)')
            snp_counter = summary_stats.shape[0]
        
        if self.remove_strand_ambig_snps:
            summary_stats['A12'] = summary_stats.A1 + summary_stats.A2
            summary_stats['A12_sorted'] = summary_stats['A12'].apply(lambda x : ''.join(sorted(x)))
            summary_stats = summary_stats[~summary_stats['A12_sorted'].isin(['AT', 'CG'])].reset_index(drop=True)
            print(f'-{snp_counter - summary_stats.shape[0]:,} of the remaining SNPs are strand-ambiguous ({summary_stats.shape[0]:,} SNPs remain)')
            snp_counter = summary_stats.shape[0]

            summary_stats['indel_status'] = summary_stats['A12'].apply(lambda x : 'N' if len(list(x)) == 2 else 'Y')
            summary_stats = summary_stats[summary_stats['indel_status'] == 'N'].reset_index(drop=True)
            print(f'-{snp_counter - summary_stats.shape[0]:,} of the remaining SNPs are indels ({summary_stats.shape[0]:,} SNPs remain)')
            snp_counter = summary_stats.shape[0]

            summary_stats.drop(columns=['A12', 'A12_sorted', 'indel_status'], inplace=True)

        imput_info_filt = summary_stats.INFO >= self.min_imputation_info_score
        summary_stats = summary_stats[imput_info_filt].reset_index(drop=True)
        print(f'-{snp_counter - summary_stats.shape[0]:,} of the remaining SNPs have imputation scores lower than {self.min_imputation_info_score}\n')
        print(f'{summary_stats.shape[0]:,} SNPs remain')
        
        self.FRAME = summary_stats.set_index('SNP')
        self.DICT = self.FRAME.to_dict(orient='index') # turn dataframe into dict

