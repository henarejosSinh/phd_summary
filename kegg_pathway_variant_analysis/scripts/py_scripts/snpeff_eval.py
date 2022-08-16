# ihc.europa@gmail.com
# dict with (effects field) impact

class snpeff_eval:
    '''
    creates a object with the impact values from snpEff fields
    self.high = "HIGH"
    self.moderate = "MODERATE"
    self.low = "LOW"
    self.modifier = "MODIFIER"
    '''
    
    def __init__(self, value):
        self.value = value
        
