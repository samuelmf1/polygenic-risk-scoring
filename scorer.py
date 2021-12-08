import math

class ScoreSNP:

    SNP_COMPLEMENTS = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    def __init__(self, base_data, target_rsid, target_alleles) -> None:
        self.target_rsid = target_rsid
        self.alleles = list(target_alleles)
        self.alleles_str = self.connect(sorted(self.alleles))

        # Verify two SNP alleles
        if len(self.alleles) != 2 and self.alleles[0] != 'N': # how do we handle indels?
            raise Exception(f'{len(self.alleles_str)} ALLELES PRESENT (INDEL?)')
        
        self.reversed = False # initialize strand reversal indicator boolean to False

        # Attempt to gather PGC Stats for the inputted SNP
        try:
            self.base_snp_dict = base_data[self.target_rsid]
            self.cyto_loc = 'chr' + str(self.base_snp_dict['CHR']) + ':' + str(self.base_snp_dict['BP'])
            self.risk_allele = self.base_snp_dict['A1']
            self.alt_allele = self.base_snp_dict['A2']
            
            """ Compares the inputted alleles to the SNP's risk allele found in the PGC data and updates alleles 
            to their complements, if necessary, to align all values. Note: strand ambiguous SNPs are invalid
            and raise an except in the init method. """
            if self.risk_allele not in self.alleles and self.risk_allele in self.__complementary_alleles__():
                self.alleles = self.__complementary_alleles__()
                self.alleles_str = self.connect(self.alleles)
                self.reversed = True

        except KeyError:
            raise KeyError('SNP NOT IN PGC DATA')

        self.calculate_score()

    def __str__(self) -> str:
        """ Returns a string representation of the class attributes. """
        components_str = [
            f'{self.target_rsid} ({self.cyto_loc})',
            f'Discov alleles:\t{self.risk_allele}{self.alt_allele} ({self.risk_allele}=risk)',
            f'Target alleles:\t{self.alleles_str} (reversed? {self.reversed})',
            f'Risk allele count: {self.risk_allele_count}',
            f'Odds ratio: {round(self.odds_ratio,3)} => ln odds ratio: {round(self.log_odds_ratio,3)}',
            f'Risk score: {round(self.risk_allele_count,3)} * {round(self.log_odds_ratio,3)} = {round(self.score,4)}'
        ]
        
        return_str = ''
        for c in components_str:
            return_str += c + '\n'
        return return_str
        
    @staticmethod
    def connect(item) -> str:
        return ''.join(item)

    def __complementary_alleles__(self) -> list:
        """ Returns a list of complementary alleles. """
        return [self.SNP_COMPLEMENTS[a] for a in self.alleles]

    def __verify_strand_alignment__(self) -> None:
        """ Compares the inputted alleles to the SNP's risk allele found in the PGC data and updates alleles 
            to their complements, if necessary, to align all values. Note: strand ambiguous SNPs are invalid
            and raise an except in the init method. """
        if self.risk_allele not in self.alleles and self.risk_allele in self.__complementary_alleles__():
            self.alleles = self.__complementary_alleles__()
            self.alleles_str = self.connect(self.alleles)
            self.reversed = True

    def calculate_score(self) -> None:
        """ Calculates risk score for the inputted allele. """
        self.odds_ratio = self.base_snp_dict['OR']
        self.log_odds_ratio = math.log(self.odds_ratio) # natural log (ln)

        """ Filters down list of self.alleles for matches with the risk allele. 
            Returns the length of the final list (i.e., number of risk alleles present). """
        self.risk_allele_count = len(list(filter(lambda a : a == self.risk_allele, self.alleles)))

        self.score = self.log_odds_ratio * self.risk_allele_count