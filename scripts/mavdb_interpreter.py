#!/usr/bin/env python
# mavedb_score_interpreter.py

class MaveDBScoreInterpreter:
    """
    A class for interpreting MaveDB scores based on different experimental assays.
    
    This class contains methods for different score interpretation paradigms and 
    handles context-dependent protein stability where some proteins function 
    better with less stability.
    """
    
    def __init__(self):
        # Dictionary mapping interpretation keywords to methods
        self.interpretation_map = {
            # K50 value patterns
            "k50": self.interpret_k50,
            "log10 k50": self.interpret_k50,
            "resistance to digestion": self.interpret_k50,
            "protease digestion": self.interpret_k50,
            "trypsin digestion": self.interpret_k50,
            "chymotrypsin digestion": self.interpret_k50,
            
            # dG value patterns
            "dg value": self.interpret_dg,
            "gibbs free energy": self.interpret_dg,
            "free energy": self.interpret_dg,
            "energy of folding": self.interpret_dg,
            "energy of unfolding": self.interpret_dg,
            
            # Log-transformed values
            "log2": self.interpret_log2,
            
            # Normalized scores around 0 for WT (not 1)
            "normalized such that wt = 0": self.interpret_normalized_zero,
            "scores are normalized such that wt = 0": self.interpret_normalized_zero,
            "scores are normalized around 0": self.interpret_normalized_zero,
            "scores are centered around zero": self.interpret_normalized_zero,
            
            # Normalized around percentiles
            "2.5th percentile": self.interpret_percentile_norm,
            
            # Electrophysiology-specific
            "half-activation voltage": self.interpret_electrophysiology,
            "half-activated voltage": self.interpret_electrophysiology,
            "v1/2": self.interpret_electrophysiology,
            
            # Specialized cases
            "yeast complementation": self.interpret_yeast_complementation,
            "binding affinity": self.interpret_binding_affinity,
            "enzyme activity": self.interpret_enzyme_activity,
            "acmg classification": self.interpret_acmg,
            
            # Default hypothesis (when no specific pattern is found)
            "default": self.interpret_default
        }
        
        # Proteins that typically rely on flexibility/instability for function
        # Higher stability could be detrimental to function in these cases
        self.flexibility_dependent_proteins = {
            "intrinsically disordered": True,
            "calmodulin": True,
            "hsp90": True,
            "p53": True,
            "myosin": True,
            "kinesin": True,
            "dynein": True,
            "ion channel": True,
            "gpcr": True,
            "g protein-coupled receptor": True,
            "transcription factor": True,
            "kinase": True,
            "phosphatase": True,
            "antibody": True,
            "immunoglobulin": True,
            "protease": True,
            "enzyme": True,  # Many enzymes need flexibility for catalysis
            "chaperone": True,
            "heat shock protein": True,
            "viral fusion protein": True,
            "allosteric": True,  # Any protein described as allosteric
        }
    
    def is_flexibility_dependent(self, description):
        """
        Determine if the protein likely relies on flexibility/instability for function
        based on keywords in the description.
        
        Args:
            description (str): The score interpretation description
            
        Returns:
            bool: True if protein likely relies on flexibility for function
        """
        if not description:
            return False
            
        description = description.lower()
        
        for keyword in self.flexibility_dependent_proteins:
            if keyword in description:
                return True
                
        return False
    
    def determine_interpretation_method(self, description):
        """
        Determine which interpretation method to use based on the score description.
        
        Args:
            description (str): The Score_Interpretation field from MaveDB
            
        Returns:
            method: The appropriate interpretation method
        """
        if not description:
            return self.interpret_default
            
        description = description.lower()
        
        # Check for matches in the interpretation map
        for keyword, method in self.interpretation_map.items():
            if keyword in description:
                return method
                
        # Default to the general hypothesis if no specific pattern is found
        return self.interpret_default
    
    def interpret_score(self, score, description, pvalue=None, high_conf=None):
        """
        Main method to interpret a MaveDB score based on its description.
        
        Args:
            score (float): The MaveDB score value
            description (str): The score interpretation text
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag ("1" or "0")
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        # Determine which interpretation method to use
        method = self.determine_interpretation_method(description)
        
        # Check if this is a protein that relies on flexibility/instability
        flexibility_dependent = self.is_flexibility_dependent(description)
        
        # Apply the appropriate interpretation method
        result = method(score, pvalue, high_conf)
        
        # Adjust interpretation if protein relies on flexibility for function
        if flexibility_dependent and 'stability' in result:
            if result['stability'] == 'increased':
                # For these proteins, increased stability may actually reduce function
                result['flexibility_dependent'] = True
                result['potential_effect'] = 'may_reduce_function'
            elif result['stability'] == 'decreased':
                # Moderate decrease might enhance function in some cases
                result['flexibility_dependent'] = True
                result['potential_effect'] = 'may_enhance_function' if score > 0.5 else 'likely_deleterious'
                
        return result
    
    def interpret_default(self, score, pvalue=None, high_conf=None):
        """
        Default interpretation based on general hypothesis:
        - Score < 1: Loss of function
        - Score > 1: Gain of function
        - Score < 0: Dominant negative
        
        Args:
            score (float): MaveDB score
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag ("1" or "0")
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'lof': False,  # Loss of function
            'gof': False,  # Gain of function
            'dn': False,   # Dominant negative
            'benign': False,  # Benign/neutral
            'mavedb_ps3': False,  # Pathogenic evidence
            'mavedb_bs3': False,  # Benign evidence
        }
        
        # Dominant negative
        if score < 0:
            result['dn'] = True
            result['lof'] = True
            
            if high_conf == "1":
                result['mavedb_ps3'] = True
        
        # Loss of function
        elif score < 0.8:
            result['lof'] = True
            
            if high_conf == "1" and score < 0.2:
                result['mavedb_ps3'] = True
            elif pvalue is not None and pvalue < 0.05:
                result['mavedb_ps3'] = True
        
        # Near wild-type
        elif 0.8 <= score <= 1.2:
            result['benign'] = True
            
            if high_conf == "1":
                result['mavedb_bs3'] = True
            elif pvalue is not None and pvalue > 0.05:
                result['mavedb_bs3'] = True
        
        # Gain of function
        else:  # score > 1.2
            result['gof'] = True
            
        return result
    
    def interpret_k50(self, score, pvalue=None, high_conf=None):
        """
        Interpret protein stability scores based on K50 values.
        Higher log10 K50 values indicate greater stability (resistance to digestion).
        
        Args:
            score (float): log10 K50 value
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'stability': None,
            'potential_effect': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        # Define thresholds for stability changes (these should be calibrated based on data)
        # For K50 values, higher scores mean greater stability
        if score < -0.5:  # Significantly less stable
            result['stability'] = 'significantly_decreased'
            result['potential_effect'] = 'likely_deleterious'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif -0.5 <= score < 0:  # Moderately less stable
            result['stability'] = 'moderately_decreased'
            result['potential_effect'] = 'possibly_deleterious'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif 0 <= score <= 0.5:  # Near wild-type stability
            result['stability'] = 'unchanged'
            result['potential_effect'] = 'likely_neutral'
            
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
            
        elif 0.5 < score <= 1:  # Moderately more stable
            result['stability'] = 'moderately_increased'
            result['potential_effect'] = 'possibly_beneficial'
            
        else:  # score > 1, Significantly more stable
            result['stability'] = 'significantly_increased'
            result['potential_effect'] = 'possibly_beneficial'
            
        return result
    
    def interpret_dg(self, score, pvalue=None, high_conf=None):
        """
        Interpret protein stability scores based on dG values.
        LOWER dG values (more negative) indicate GREATER stability.
        This is the opposite relationship from K50 values.
        
        Args:
            score (float): dG value or related measure
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'stability': None,
            'potential_effect': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        # For dG values, LOWER scores mean GREATER stability (opposite of K50)
        if score > 0.5:  # Significantly less stable
            result['stability'] = 'significantly_decreased'
            result['potential_effect'] = 'likely_deleterious'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif 0 < score <= 0.5:  # Moderately less stable
            result['stability'] = 'moderately_decreased'
            result['potential_effect'] = 'possibly_deleterious'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif -0.5 <= score <= 0:  # Near wild-type stability
            result['stability'] = 'unchanged'
            result['potential_effect'] = 'likely_neutral'
            
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
        
        elif -1 <= score < -0.5:  # Moderately more stable
            result['stability'] = 'moderately_increased'
            result['potential_effect'] = 'possibly_beneficial'
            
        else:  # score < -1, Significantly more stable
            result['stability'] = 'significantly_increased'
            result['potential_effect'] = 'possibly_beneficial'
            
        return result
    
    def interpret_log2(self, score, pvalue=None, high_conf=None):
        """
        Interpret log2-transformed values:
        - For log2 scores, 0 often represents wild-type
        - Negative values represent decreased effect
        - Positive values represent increased effect
        
        Args:
            score (float): log2 score
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'effect': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        # Log2 scores typically center around 0 for wild-type
        if score < -1:
            result['effect'] = 'strongly_decreased'
            
            # With high confidence or significant p-value
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif -1 <= score < -0.5:
            result['effect'] = 'moderately_decreased'
            
            # With high confidence or significant p-value
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif -0.5 <= score <= 0.5:
            result['effect'] = 'unchanged'
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
            
        elif 0.5 < score <= 1:
            result['effect'] = 'moderately_increased'
            
        elif score > 1:
            result['effect'] = 'strongly_increased'
            
        return result
    
    def interpret_normalized_zero(self, score, pvalue=None, high_conf=None):
        """
        Interpret scores normalized such that WT = 0 (not 1).
        Often the 2.5th percentile = -1.
        Positive values indicate better performance than WT.
        
        Args:
            score (float): Normalized score
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'effect': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        # WT = 0, more negative = worse than WT
        if score < -1:
            result['effect'] = 'strongly_decreased'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif -1 <= score < -0.5:
            result['effect'] = 'moderately_decreased'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif -0.5 <= score <= 0.5:  # Near WT
            result['effect'] = 'unchanged'
            
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
            
        elif 0.5 < score <= 1:
            result['effect'] = 'moderately_increased'
            
        elif score > 1:
            result['effect'] = 'strongly_increased'
            
        return result
    
    def interpret_percentile_norm(self, score, pvalue=None, high_conf=None):
        """
        Interpret scores normalized where WT = 0 and 2.5th percentile = -1.
        These scores are normalized relative to the distribution of all variants.
        
        Args:
            score (float): Percentile-normalized score
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'effect': None,
            'percentile_interpretation': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        if score < -1:
            result['effect'] = 'strongly_decreased'
            result['percentile_interpretation'] = 'below_2.5th_percentile'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif -1 <= score < -0.5:
            result['effect'] = 'moderately_decreased'
            result['percentile_interpretation'] = 'between_2.5th_and_16th_percentile'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif -0.5 <= score <= 0.5:
            result['effect'] = 'unchanged'
            result['percentile_interpretation'] = 'near_median'
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
            
        elif 0.5 < score <= 1:
            result['effect'] = 'moderately_increased'
            result['percentile_interpretation'] = 'between_84th_and_97.5th_percentile'
            
        elif score > 1:
            result['effect'] = 'strongly_increased'
            result['percentile_interpretation'] = 'above_97.5th_percentile'
            
        return result
    
    def interpret_electrophysiology(self, score, pvalue=None, high_conf=None):
        """
        Interpret electrophysiology scores like half-activation voltage (V1/2).
        
        Args:
            score (float): Electrophysiology measurement (e.g., V1/2)
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'channel_property': None,
            'potential_effect': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        # For half-activation voltage, interpretation depends on channel type
        # More negative: Channel activates at more hyperpolarized potentials
        # More positive: Channel activates at more depolarized potentials
        if score < -0.5:
            result['channel_property'] = 'activates_at_more_negative_voltages'
            result['potential_effect'] = 'possibly_gain_of_function'
            
        elif -0.5 <= score <= 0.5:
            result['channel_property'] = 'unchanged_activation'
            result['potential_effect'] = 'likely_neutral'
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
            
        else:  # score > 0.5
            result['channel_property'] = 'activates_at_more_positive_voltages'
            result['potential_effect'] = 'possibly_loss_of_function'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        return result
    
    def interpret_yeast_complementation(self, score, pvalue=None, high_conf=None):
        """
        Interpret yeast complementation assays where higher scores indicate
        better complementation of the yeast phenotype.
        
        Args:
            score (float): Yeast complementation score
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'complementation': None,
            'potential_effect': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        # For yeast complementation, often:
        # 0 = null controls, 1 = wild-type controls
        if score < 0.2:
            result['complementation'] = 'severely_impaired'
            result['potential_effect'] = 'likely_loss_of_function'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif 0.2 <= score < 0.8:
            result['complementation'] = 'partially_impaired'
            result['potential_effect'] = 'partial_loss_of_function'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif 0.8 <= score <= 1.2:
            result['complementation'] = 'wild_type_like'
            result['potential_effect'] = 'likely_neutral'
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
            
        else:  # score > 1.2
            result['complementation'] = 'enhanced'
            result['potential_effect'] = 'possible_gain_of_function'
            
        return result
    
    def interpret_binding_affinity(self, score, pvalue=None, high_conf=None):
        """
        Interpret binding affinity scores where higher scores typically
        indicate stronger binding.
        
        Args:
            score (float): Binding affinity score
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'binding': None,
            'potential_effect': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        # For binding assays, interpretation depends on context
        if score < 0.5:
            result['binding'] = 'significantly_reduced'
            result['potential_effect'] = 'likely_loss_of_function'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif 0.5 <= score < 0.8:
            result['binding'] = 'moderately_reduced'
            result['potential_effect'] = 'partial_loss_of_function'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif 0.8 <= score <= 1.2:
            result['binding'] = 'wild_type_like'
            result['potential_effect'] = 'likely_neutral'
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
            
        elif 1.2 < score <= 1.5:
            result['binding'] = 'moderately_enhanced'
            result['potential_effect'] = 'possible_gain_of_function'
            
        else:  # score > 1.5
            result['binding'] = 'significantly_enhanced'
            result['potential_effect'] = 'likely_gain_of_function'
            
        return result
    
    def interpret_enzyme_activity(self, score, pvalue=None, high_conf=None):
        """
        Interpret enzyme activity scores where higher scores typically
        indicate increased activity.
        
        Args:
            score (float): Enzyme activity score
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'activity': None,
            'potential_effect': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
        
        # For enzyme activity assays
        if score < 0.2:
            result['activity'] = 'severely_reduced'
            result['potential_effect'] = 'likely_loss_of_function'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif 0.2 <= score < 0.8:
            result['activity'] = 'moderately_reduced'
            result['potential_effect'] = 'partial_loss_of_function'
            
            if high_conf == "1" or (pvalue is not None and pvalue < 0.05):
                result['mavedb_ps3'] = True
                
        elif 0.8 <= score <= 1.2:
            result['activity'] = 'wild_type_like'
            result['potential_effect'] = 'likely_neutral'
            if high_conf == "1" or (pvalue is not None and pvalue > 0.05):
                result['mavedb_bs3'] = True
            
        elif 1.2 < score <= 2:
            result['activity'] = 'moderately_enhanced'
            result['potential_effect'] = 'possible_gain_of_function'
            
        else:  # score > 2
            result['activity'] = 'significantly_enhanced'
            result['potential_effect'] = 'likely_gain_of_function'
            
        return result
    
    def interpret_acmg(self, score, pvalue=None, high_conf=None):
        """
        Interpret scores based on ACMG classification system.
        
        Args:
            score (float): ACMG-based score
            pvalue (float, optional): P-value if available
            high_conf (str, optional): High confidence flag
            
        Returns:
            dict: Dictionary containing interpretation results
        """
        result = {
            'acmg_class': None,
            'mavedb_ps3': False,
            'mavedb_bs3': False
        }
            
        return result


def map_interpretation_to_method(interpretation_text):
    """
    Map an interpretation text to the appropriate method in the MaveDBScoreInterpreter class.
    
    Args:
        interpretation_text (str): The score interpretation description
        
    Returns:
        str: The name of the most appropriate interpretation method
    """
    interpreter = MaveDBScoreInterpreter()
    method = interpreter.determine_interpretation_method(interpretation_text)
    method_name = method.__name__
    
    # Check if this is a flexibility-dependent protein
    flexibility_dependent = interpreter.is_flexibility_dependent(interpretation_text)
    if flexibility_dependent:
        return f"{method_name} (context: flexibility-dependent protein)"
    
    return method_name


def get_score_classification(score_value, interpretation_text, pvalue=None, high_conf=None):
    """
    Helper function to classify a MaveDB score using the interpreter class.
    
    Args:
        score_value (float): The numerical score from MaveDB
        interpretation_text (str): The score interpretation description
        pvalue (float, optional): P-value if available
        high_conf (str, optional): High confidence flag
        
    Returns:
        dict: Dictionary containing classification results
    """
    interpreter = MaveDBScoreInterpreter()
    return interpreter.interpret_score(score_value, interpretation_text, pvalue, high_conf)