# ihc.europa@gmail.com
# class to calculate and return effect score 


# import snpeff_eval


class EffectScore:
    '''
    Implements a class holding the necessary variables to calculate 
    the damage of certain variants in a node:
    string effect
    snpeff_eval object impact
    float score
    int hash
    '''
    
    def __init__(self, effect, impact, score):
        self.effect = effect
        self.impact = impact
        self.score = score
        self.dict = {}

    @staticmethod
    def __hash__(self):
        return(hash(self))
        
    @staticmethod
    def __create_obj(effect, impact, score):
        return EffectScore(effect, impact, score)
    
    def hash_dict(self):
        hash_val = self.__hash__()
        self.dict = {hash_val: self.effect}
        
    # get methods
    def get_effect(self):
        return self.effect
    
    def get_impact(self):
        return self.impact
    
    def get_score(self):
        return self.score
    
    def get_dict(self):
        return self.dict
    
    def to_string(self):
        return f'Effect Score: effect= {self.effect}, impact= {self.impact}score= {self.score}'

    # compare methods
    # compare if get_score() and effect_score.get_score() are equal
    # compare if self.effect == effect_score().effect