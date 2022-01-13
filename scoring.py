import math
from tqdm import tqdm

class ScoreEngine:

    SNP_COMPLEMENTS = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    def __init__(self, base_data) -> None:
        if isinstance(base_data, dict):
            self.base_data = base_data
        else:
            self.base_data = base_data.to_dict(orient='index')

    def score_polygenic_risk(self, target_dict:dict) -> float:
        
        target = list(target_dict.items())
        prs_cum = 0

        for index in tqdm(range(len(target_dict))):
            rsid, alleles = target[index]
            prs_cum += self.score_snp(rsid, alleles)

        return prs_cum

    def score_snp(self, target_rsid:str, target_alleles:str, **kwargs) -> float:

        alleles = list(target_alleles)
        alleles_str = ''.join(sorted(alleles))
        reversed = False # initialize strand reversal indicator boolean to False

        # Verify two SNP alleles
        if len(alleles) != 2 and alleles[0] != 'N': # how do we handle indels?
            raise Exception(f'{len(alleles_str)} ALLELES PRESENT (INDEL?)')
        
        # Attempt to gather PGC Stats for the inputted SNP
        try:
            base_snp_dict = self.base_data[target_rsid]
            cyto_loc = 'chr' + str(base_snp_dict['CHR']) + ':' + str(base_snp_dict['BP'])
            risk_allele = base_snp_dict['A1']
            alt_allele = base_snp_dict['A2']
            
            """ Compares the inputted alleles to the SNP's risk allele found in the base data and updates alleles 
            to their complements, if necessary, to align all values. Note: strand ambiguous SNPs are invalid
            and raise an except in the init method. """
            complementary_alleles = [self.SNP_COMPLEMENTS[a] for a in alleles]
            if risk_allele not in alleles and risk_allele in complementary_alleles:
                alleles = complementary_alleles
                alleles_str = ''.join(alleles)
                reversed = True
        except KeyError:
            raise('SNP not in base data')

        """ Calculates risk score for the inputted allele. """
        odds_ratio = base_snp_dict['OR']
        log_odds_ratio = math.log(odds_ratio) # natural log (ln)

        """ Filters down list of self.alleles for matches with the risk allele. 
            Returns the length of the final list (i.e., number of risk alleles present). """
        risk_allele_count = len(list(filter(lambda a : a == risk_allele, alleles)))

        score = log_odds_ratio * risk_allele_count

        if 'print' in kwargs and kwargs['print'] is True:
            components_str = [
                f'{target_rsid} ({cyto_loc})',
                f'Discov alleles:\t{risk_allele}{alt_allele} ({risk_allele}=risk)',
                f'Target alleles:\t{alleles_str} (reversed? {reversed})',
                f'Risk allele count: {risk_allele_count}',
                f'Odds ratio: {round(odds_ratio,3)} => ln odds ratio: {round(log_odds_ratio,3)}',
                f'Risk score: {round(risk_allele_count,3)} * {round(log_odds_ratio,3)} = {round(score,4)}'
            ]
            return_str = ''
            for c in components_str:
                return_str += c + '\n'
            print(return_str)

        return score