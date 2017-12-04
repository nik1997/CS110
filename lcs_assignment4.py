import numpy as np
import pandas as pd

def longest_subsequence(string1, string2):
    '''
        Returns the length of longest subsequence between two strings.
        
        Input:
            string1, string2 = Two strings whos LCS is to be found.
            
        Output:
            The length of longest subsequence between two strings
    '''
    
    table = np.zeros([len(string1)+1, len(string2)+1])
        # Stores the length of LCS till each first few elements
        # +1 done to add one extra column and row to have first row and column as empty
     
    for i in range(len(string1) + 1):
        # Selecting each character of first string and comparing it with each character of second string
        for j in range(len(string2) + 1):
            if i == 0 or j == 0: 
                table[i,j] = 0
            elif string1[i-1] == string2[j-1]: 
                # If two selected elements are same, adding the length of LCS till that point 
                table[i,j] = table[i-1, j-1] + 1
            else:
                # If they are not same, updating length of LCS of at point with the last maximum value.
                table[i,j] = max(table[i, j-1], table[i-1, j])
    
    return table[len(string1), len(string2)]
    
def length_table_creator(genes):
    '''
        Creates the table with LCS for each pair of strings.
        
        Input:
            genes = A 1-d array with a set in each element. The 1st element of set is its index and 2nd element is the gene-string.
        
        Output:
            len_LCS_table = A 2-d array that stores length of LCS of each pair. 
        
    '''
    len_LCS_table = np.zeros([len(genes), len(genes)])  
    for i in genes:
        for j in genes:
            if i == j:
                continue
            else:
                #Check if the LCS of a certain pair of string have found or not. Only calling the longest_subsequence funtion if required.
                if len_LCS_table[j[0],i[0]] == 0:
                    len_LCS_table[i[0],j[0]] = longest_subsequence(i[1], j[1])
                else:
                    len_LCS_table[i[0], j[0]] = len_LCS_table[j[0], i[0]]
    
    return len_LCS_table



'''
    Base Rules and Assumptions:
    - Each gene can only have two offsprings
    - If LCS of the difference between length of LCS of one gene with another genes is less than or equals to a certain relationship factor(default used as 3) and the length of maximum LCS, the gene has
    direct relationship with those genes
        - For instance, if Str 1 has LCS of length 83 and 81 with Str 4, it has a direct relationship with both of them.
        - Direct relationship means either str 1 is children or parent of str 4
    - If gene has only two LCS length within the difference of 3 from the maximum LCS, the gene is the base gene
    - If gene has no LCS with difference of length between maximum LCS gene and LCS with other genes, the gene is the last child gene.
    
    Glossary:
    direct_relationship: means that two genes have a parent-child relationship
        A has direct_relationship with B, if the difference between maximum length of LCS between A and any other gene
        and, the length of LCS between A and B is less than or equals to 5.
        This has been defined, based on the analysis of the pattern in lengths of LCS
            
'''
class rltn_tree:
    '''
        @params:
            len_LCS_table = A 2-d array that stores length of LCS of each pair of strings
    '''
    
    def __init__(self, len_LCS_table):
        self.len_LCS_table = len_LCS_table
        self.rltn_factor = 5 # The value by which one length of LCS should be smaller than max LCS to have direct relationship
        self.ordered_array = [] #Stores genes in order of their generation
        ordered_array = self.rltn_tree_maker() # Creating the ordered_array with relationship when the class is used.
        
    def find_children(self, parent, gene):
        '''
            Inputs:
                parent : The parent of the gene whose children are to be found.
                gene : The gene whose children are to be found

            Output:
                No output but it updates the ordered_array that stores the genes in order of their generation

            Process:
                Finds the children genes who have direct relationship with the gene
                If there are two children genes, then prints the children genes
                Then, calls the find_children function again by passing the children gene
        '''

        max_LCS_length = max(self.len_LCS_table[gene, :]) #Finding maximum length of gene

        genes_children = [] # Array to store the genes id of the children

        for j in range(len(self.len_LCS_table[gene])):
            # Two genes have direct relationship if follows the following condition
            if max_LCS_length - self.len_LCS_table[gene, j] <= self.rltn_factor:
                # If the other gene is not parent, then the other gene should be the children
                if j != parent:
                    genes_children.append(j)
        if len(genes_children) == 2:
            # Adding the children genes to the ordered_array
            self.ordered_array.append(genes_children[0])
            self.ordered_array.append(genes_children[1])

            # Calling the recursive function to find the children of children
            self.find_children(gene, genes_children[0])
            self.find_children(gene, genes_children[1])
        else:
            # If the gene has only one other gene with direct relation, then it will end the recursion and return
            return

    def rltn_tree_maker(self):
        '''
            Inputs:
                len_LCS_table : A 2x2 array that has the length of LCS of all pairs of genes

            Outputs:
                ordered_array : A 1-D array that stores elements in their order of generation 

            Process:
                - First finds the base gene by finding the number of genes it has direct relationship with.
                - If it has direct relationship with only 2, it is the base gene.
                - Then, calls find_children function to find the children of the children of the base genes           
        '''
        
        self.ordered_array = [] #Initializing ordered_array to [] to make sure it empty before another ordered tree is made 
        
        for i in range(len(self.len_LCS_table)): #Checking all genes available
            max_LCS_length = max(self.len_LCS_table[i,:])
            genes_children = [] #Stores the children genes

            for j in range(len(self.len_LCS_table[i,:])):
                if max_LCS_length - self.len_LCS_table[i, j] <= self.rltn_factor:
                    genes_children.append(j)

            if len(genes_children) == 2:
                # Breaks the loop if there are only 2 children genes, since this is the base gene
                break

        # Appending the Base Parent's Gene and its Children's gene
        self.ordered_array.append(i)
        self.ordered_array.append(genes_children[0])
        self.ordered_array.append(genes_children[1])

        # Calling to find childrens of each child
        self.find_children(i, genes_children[0])
        self.find_children(i, genes_children[1])
        
        return self.ordered_array

    def rltn_tree_printer(self):
        '''
            Process:
                Prints each generation in a tree order.
                Since we assume that one gene can only make two offsprings, 
                we increase the last index and loop by multiplying it by two, so that it can accomodate the
                exponential growth

            Input:
                ordered_array = A 1-d array with genes in order of their generation

            Output:
                Prints each generation in a separate line 
        '''
        try:
            if self.ordered_array == []:
                raise Exception
            start_index = 0
            last_index = 1
            loop = True
            while loop:
                if last_index <= len(self.ordered_array):
                    for i in range(start_index,last_index):
                        print self.ordered_array[i],
                else:
                    loop = False
                    for i in range(start_index,len(self.ordered_array)):
                        print self.ordered_array[i],
                print "\n"
                start_index = last_index
                last_index = last_index * 2 + 1   
        except:
            print "Ordered_Array is empty. Please run rltn_tree_maker before running this function"
            
            
# Testing 
genes = [(0,'TTCTACGGGGGGAGACCTTTACGAATCACACCGGTCTTCTTTGTTCTAGCCGCTCTTTTTCATCAGTTGCAGCTAGTGCATAATTGCTCACAAACGTATC'), 
            (1,'TCTACGGGGGGCGTCATTACGGAATCCACACAGGTCGTTATGTTCATCTGTCTCTTTTCACAGTTGCGGCTTGTGCATAATGCTCACGAACGTATC'), 
            (2,'TCTACGGGGGGCGTCTATTACGTCGCCAACAGGTCGTATGTTCATTGTCATCATTTTCATAGTTGCGGCCTGTGCGTGCTTACGAACGTATTCC'), 
            (3,'TCCTAACGGGTAGTGTCATACGGAATCGACACGAGGTCGTATCTTCAATTGTCTCTTCACAGTTGCGGCTGTCCATAAACGCGTCCCGAACGTTATG'), 
            (4,'TATCAGTAGGGCATACTTGTACGACATTCCCCGGATAGCCACTTTTTTCCTACCCGTCTCTTTTTCTGACCCGTTCCAGCTGATAAGTCTGATGACTC'), 
            (5,'TAATCTATAGCATACTTTACGAACTACCCCGGTCCACGTTTTTCCTCGTCTTCTTTCGCTCGATAGCCATGGTAACTTCTACAAAGTTC'), 
            (6,'TATCATAGGGCATACTTTTACGAACTCCCCGGTGCACTTTTTTCCTACCGCTCTTTTTCGACTCGTTGCAGCCATGATAACTGCTACAAACTTC')]


def test():
    len_LCS_table = length_table_creator(genes)
    headers = [i[0] for i in genes]
    df = pd.DataFrame(len_LCS_table, index = headers, columns=headers)
    print df, "\n"
    
    print "Printing Relationship Tree \n"
    tree = rltn_tree(len_LCS_table)
    tree.rltn_tree_maker()
    tree.rltn_tree_printer()

test()
