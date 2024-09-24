import sys
class OrfFinder():
    """
    Class to represent a sequence of DNA and find gene candidates. This code has methods to find candidates on all 3 frames in the positive and negative strand and to print the results in an approriate format
    """
    
    pairs = {"A":"T","C":"G","T":"A","G":"C"}

    def __init__ (self,seq,longGene=True,minGene=300,start=["ATG"],stop=["TAA"]):
        """
        intializes a class with attributes for all extra credit options, and sequence and uses these values to find open reading frames and add them to an attribute 
        """
        self.start = start
        self.stop = stop
        self.frames = []   
        self.seq = seq
        self.minGene = minGene
        self.longGene = longGene
        self.addseq(seq,longGene,minGene,start,stop)
    
    def frameFinder(self,seq,strandpos):
        """
        Returns a list of Open reading frames on all 3 possible frames from a sequence of DNA either the positive or negative strand
        """
        start1 = self.start
        stop = self.stop
        frames = []
        for frame in range(0,3): #looping through 3 possible frames
            stopfound = False #tracks whether a stop codon has been encountered in this frame
            startpos=[] #list of start codons found
            for index in range(frame,len(seq)-2,3): #loops through 3 bases at a time using the frame# from the parent loop as an offset to search the right frame
                if seq[index:index+3] in start1: #checking for start codons 
                    if index != 0: #preventing duplicate frames when we handle hanging cases later
                        startpos.append(index) 
                if seq[index:index+3] in stop: #checking for stop coodns
                    startpos = sorted(startpos)
                    for start in startpos: #adding frames for all start positions found before this stop codon
                        if ((index+3)-start) >= self.minGene:
                            if strandpos == True:
                                frames.append([frame+1,start+1,index+3,(index+3)-start])
                            elif strandpos == False:
                                frames.append([-(frame+1),len(seq)-(index+2),len(seq)-start,(index+3)-start])
                        if self.longGene == True: #if longest gene is true we shoudl only append the frame found for the earliest start position before the stop coodn 
                            break 
                        #exit start positons loop all start positions prior to this stop handled, hanging cases adressed in block below 
                    if stopfound == False: #block to handle hanging cases with no start position, if a stop has been found before this one it can not be a no start hang because a stop cant be read through
                        if startpos == []:
                            if (index+3)>=self.minGene:
                                if strandpos == True:
                                    frames.append([(frame+1),1,index+3,index+3])
                                elif strandpos == False:
                                    frames.append([-(frame+1),index+3,len(seq),index+3])
                    #reset start positions and update stop codon tracker to found
                    startpos = []
                    stopfound=True
                    #exit stop found conditional block, all start positons prior to this stop codon and any relevant hanging cases handled
                #exit loop iterating through codons in frame, all codons in frame searched 
            if stopfound == False: #block to handle hanging cases for a frame with no stops 
                for start in startpos:
                    if strandpos == True:
                        frames.append([frame+1,start+1,len(seq),len(seq)-start,])
                    elif strandpos == False:
                        frames.append([-(frame+1),1,len(seq)-start+1,len(seq)-start])
                if startpos == []:
                        if strandpos == True:
                            frames.append([(frame+1),1,len(seq),len(seq)])
                        elif strandpos == False:
                            frames.append([-(frame+1),1,len(seq),len(seq)])
            #exit parent loop, iterate to next frame if relevant  
        frames = sorted(frames, key=lambda x:(-x[3])) 
        return frames
    
    def findReverse(self,str):
        """
        Returns the complimentary seqeunce to a strand in the 5'-3' direction

        input: ACTG
        Ouput: CAGT
        """
        compliment = ""
        for x in str:
            for char in x:
                compliment += self.pairs[char]
        return compliment[::-1]
    
    def addseq(self,seq,longGene,minGene,start,stop):
        """
        Takes in a sequence and adds the frames on both reverse and forward strand to the frames attributes
        """
        seqReverse = self.findReverse(seq) #find reverse compliment in 5-3 directio
        self.frames = self.frames + (self.frameFinder(seq,True)) #find frames on positive strand
        self.frames = self.frames + (self.frameFinder(seqReverse,False)) #find frames on reverse strand 
    
    def printer(self):
        """
        prints all frames in appropriate format 
        """

        frames = sorted(self.frames, key=lambda x: -x[3])
        
        max_length = max(len(str(item[2])) for item in frames)

        # Print the items with consistent spacing
        for item in frames:
            print(f"{item[0]:+} {item[1]:<6}..{item[2]:<6}  {item[3]:<{max_length}}")
        

class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname,'r',encoding='utf-8')
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()
            
        yield header,sequence
class NucParams:
    """
    class to represent a genome with methods to return: sequence length,a list of codons, codon composition, codon relative frequency, 
    a list of amino acids coded for, a list of amino acid composition, and GC content by percentage
    """
    
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}
    
    validAA = "-ACDEFGHIKLMNPQRSTVWY"
    
    def __init__ (self, inString=''):
        '''
        Instantiate object attributes: amino acid and codon composition, nucleotide composition, and total nucleotide count all as empty. Use addseq to process and 
        add instring parameter to these 
        '''
        #self.aaComp = {"F":0,"L":0,"I":0,"M":0,"V":0,"S":0,"P":0,"T":0,"A":0,"Y":0,"-":0,"H":0,"Q":0,"N":0,"K":0,"D":0,"E":0,"C":0,"W":0,"R":0,"S":0,"G":0}
        self.aaComp = {'-':0,'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}
        self.codonComp = {codon: 0 for codon in self.rnaCodonTable.keys()}
        self.nucComp = {"A":0,"C":0,"T":0,"G":0,"U":0,"N":0}
        self.nucCount = 0
        self.addSeq(inString)
    def codonSort(self,inSeq):
        """
        Sort an input string into a list of rna codons 
        
        input: acggcu
        output ["acg","gcu"]
        """
        codonList = []
        count = 0
        codon = ""
        #for loop splits in string every 3 bases, uppercases input, and changes t to u 
        for index, baseLow in enumerate(inSeq):
            base = baseLow.upper()
            if base == "T":
                base = "U"
            if count<2:  
                codon += base
                count += 1
            else: #if the count is 2 once we finish adding this base the codon will have 3 bases, 0,1, and 2
                codon+=base
                codonList.append(codon) #codon now contains 3 bases in order, it is added to the list of codons
                #count and codon are reset
                count = 0 
                codon = ""
        #for loop searches list for invalid codons and deletes them
        #its important to use some version of this code to remove nonsense characters and codons without inducing a frameshift 
        for codon in codonList: 
            if codon not in self.rnaCodonTable:
                codonList.remove(codon)
        
        return codonList
    def addSeq (self, inSeq):
        """
        Processes sequence to find codon and amino acid composition and nucleotide composition and count and add them to the object attributes
        """
        
        #local variable instantiation block, creates local instances of a list of codons and amino acids and dictionaries of their counts for the in seqeunce
        #Instantiating local variables for input parameters and important variables taken from a class method or attribute is a personal style choice
        #I think its worth the extra line for the improved style, although its often unecessary it comes in handy enough to be worth getting in the habit of doing for me 
        
        #list method block
        #calling codonsort and passing it the input seqeunce parameter, returns a codon list
        localCodonList = self.codonSort(inSeq) 
        #calling amino acid list method and passing it the codon list from last line, returns an amino acid list
        localaaList = self.aaList(localCodonList) 
            #list methods are only necessary for local use, kept around for personal preference to keeping code in separate methods
        
        
        #dictionary method block
        #calling codon composition method and passing it the codon list, returns codon comp dictionary
        localCodonComp = self.codonComposition(localCodonList) 
        #calling amino acid composition method and passing it amino acid list from a few lines back, returns amino acid comp dictionary 
        localAAComp = self.aaComposition(localaaList) #
            #composition methods are the real meat of this method, these rely on the list methods for their inputs, improves code clarity imo
            #comes at the cost of making a bit of a labrynthine program, for future reference pay close attention to code order debug met a lot
            #of issues relating to the order the interpreter parses these methods, + a mistake in one method can screw 4 others up 
        
        #using for loops to iterate through codon, nucleotide and amino acid dictionaries and add their counts to the larger object dictionary
        for codon, value in localCodonComp.items(): 
            self.codonComp[codon]+= value #adding count of each codon in the current sequence to the count registered in the object attribute
            
        for aa, value in localAAComp.items():
            self.aaComp[aa]+= value #adding count of each amino acid in the current sequence to the count registered in the object attribute
        
        for nuc, count in self.nucComposition(inSeq).items():
            self.nucComp[nuc] += count #adding count of each nucleotide in the current sequence to the count registered in the object attribute
            
        self.nucCount += self.nucCount1(self.nucComposition(inSeq)) #adding nucleotide count of this sequence to the total 
    def aaList(self,codonList):
        """
        Returns list of amino acids from a list of codons
        
        input: ["UUU","UUA"]
        output: ["F","L"]
        """
        localCodonList = codonList #instantiating local variable for input 
        aaList = []
        for codon in localCodonList:
            aaList.append(self.rnaCodonTable[codon])
        return aaList
    def aaComposition(self,inAAlist):
        """
        Return dictionary of amino acid counts from a list of amino acids 
        
        input: ["F","L"]
        output: {"F":1,"L":1}
        """
        protein = inAAlist #instantiating local variable for input 
        aaComp = {} 
        for aminoAcid in self.validAA: 
            if aminoAcid in protein:
                aaComp[aminoAcid]=protein.count(aminoAcid)
        return aaComp
    def nucComposition(self,inSeq):
        """
        Return dictionary of nucleotide count from a dna/rna sequence 
        input: "AGU"
        output: {"A":1,"G":1,"U":1}
        """
        #instantiating local variable for input 
        localSeq = inSeq
        nucComp = {}  
        for base in localSeq:
            nucComp[base] = localSeq.count(base)
        return nucComp
    def codonComposition(self,codonList):
        """
        Return a dictionary of codon counts from list of codons
        
        input: ["AGU","AUG"]
        Return: {"AGU":1,"AUG":1}
        """
        #instantiating local variable for input 
        localCodonList = codonList
        codonComp={}
        for item in localCodonList:
            codonComp[item]=localCodonList.count(item)
        return codonComp
    def codonFreq(self):
        """
        Return Dictionary of relative codon frequencies 
        
        input: ['GCU','GCC']
        Output: {'GCU':50,'GCC':50} 
        """
        codonFreq = {}
        #instantiating local variable 
        codonCount = self.codonComp
        aaCount = self.aaComp
        for codon, aa in self.rnaCodonTable.items():
            if codonCount[codon] != 0:
                codonFreq[codon] = (codonCount[codon]/aaCount[aa])*100
        return codonFreq
    def nucCount1(self,inDict):
        """
        Return total sequence length inncluding n bases
        
        input: "ACGT"
        output: 4
        """
        #instantiating local variable for input 
        localNucComp = inDict 
        return sum(localNucComp.values())
    def gcContent(self):
        """
        Return GC content of sequence
        
        input:"ATGC"
        Output: 50
        """
        return ((self.nucComp["G"]+self.nucComp["C"])/self.nucCount)*100
    def printer(self):
        """
        Return print statement as string containing the sequence length in Mb, the Gc content, and the relative codon bias of a sequence 
        """
        
        #instantiating local variable for input 
        codonFreq = self.codonFreq()
        outStr = ""
        #for loop iterates through alphabetically ordered amino acids checking for the codon count and freqeuncies assoicated with codons that code for
        #this acid, this task necessitates the clunky nested loop to find all associated codons, values are added to an outstring which is returned
        for aa in self.validAA:  
            for key, value in self.rnaCodonTable.items():
                if value == aa:
                    if self.codonComp[key] != 0:
                        #outStr += "\n" + str(key) + " : " + str(value) + " " + str(round(codonFreq[key],2)) + " ( " + str(self.codonComp[key]) + " )"
                        outStr += "\n" + ('{:s}:{:s} {:5.1f} ({:6d})'.format(key,value,codonFreq[key],self.codonComp[key]))
        return ("Sequence length = " + str(self.nucCount/1000000) + "mb" + "\n " + "\nGC Content = " + str(self.gcContent()) + outStr)
class ProteinParam(str):
    """Class to represent a protein"""

# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
            'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
            'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
            'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
            'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
            }

    aa = 'ACDEFGHIKLMNPQRSTVWY'
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        """
        Intialize a protein object with a dictionary that holds the count of each valid amino acid and a raw string variable that contains the input
        string washed of nonsense characters and corrected to be uppercase.
        """

        #initializing raw protein string variable to make molecular weight calculation easier
        washedProteinString = "" #pre instantiating return value

        for aa in protein: #for loop iterates through string comparing values to built in amino acid dictionary
            if aa.upper() in self.aa2mw:
                washedProteinString += aa.upper() #values that meet the conditional are uppercased and added to the raw string

        self.washedProteinString = washedProteinString

        #dictionary with amino acid count built in
        self.protein = self.aaComposition()

    def aaCount (self):
        """
        Return the total number of amino acids in the string.

        input:VLSPADKTNVKAAW
        output:14
        """
        sumAA = 0 #pre instantiating return value

        for key, value in self.protein.items(): #for loop iterates through dictionary of amino acid counts
            sumAA += value #count of each amino acids is summed to return the total number of amino acids

        return sumAA

    def pI (self,precision=.01):
        """
        Return the ph at which the net charge is the closest to 0, this is the theoretical isoelectric point.

        input:VLSPADKTNVKAAW
        output:9.88
        """
        #using a binary search method to return the theoretical isoelectric point to a specified precision. The loop
        #iterates through the range of possible ph values by searching the center point of increasingly smaller ranges.
        #if the current count(the center point) is less than 0 we know based off the shape of a titration curve that
        #the ph must be too high and any value higher does not need to be searched so the upper bound is set to the current
        #range and the loop iterates again, vice versa for values greater than 0. The loop exits when the size of the range
        #between the upper and lower bound is less than the current count +/- specified precision

        lowPH = 0
        highPH = 14
        #print(f"Precision is {precision}")
        while highPH - lowPH >= 2*precision:
            count = (highPH + lowPH)/2
            #print("count is: " + str(count))
            if self._charge(count) < 0:
                #print(f"The charge is less than 0, it is {self._charge(count)}")
                highPH = count

            elif self._charge(count) > 0:
                #print(f"The charge is greater than 0, it is {self._charge(count)}")
                lowPH = count
            else:
                #print("count")
                break
        return count

    def aaComposition (self) :
        """
        Return a dictionary where keys are the amino acids present in the string and values are the number of each amino acid present in the string.

        input:AAGG
        return:{A:2,C:0,D:0...G:2...Y:0}
        """
        protein = self.washedProteinString #local variable initalized to the value of object attirbute for readability
        aaCountDict = {} #pre instantiating return value

        for aminoAcid in self.aa2mw: #for loop iterates through string comparing values to the built in amino acid dictionary
            if aminoAcid in protein:
                aaCountDict[aminoAcid]=protein.count(aminoAcid) #values that meet the condition are added as a key, their corresponding value is their count in the string
            else:
                aaCountDict[aminoAcid]=0 #values that are not present in the string are added as a key with their corresponding value being 0
        return aaCountDict

    def _charge (self,pH):
        """
        Return the charge of the protein at a given ph value by summing the positive and negative charges of the charged residues multiplied by their count
        and adding the positive and negative charge of the c and n terminus.

        input: VLSPADKTNVKAAW, 7
        outpit 0.9980759553300174
        """
        posChargeTotal = 0 #local variable pre instantiation
        negChargeTotal = 0 #local variable pre instantiation

        #we can use a loop to check for positive and negative charged residues and use the dictionaries for the pka of both respectively to calculate total positive and negative charge
        for key, value in self.protein.items(): #iterating through amino acid count dictionary
            if key in self.aa2chargePos: #math for positively charged amino acid residues starts here
                posChargeTotal += value*(10**self.aa2chargePos[key])/((10**self.aa2chargePos[key])+(10**pH)) #using equation provided to return charge from a given residue * its count
            if key in self.aa2chargeNeg: #math for negatively charged amino acid residues starts here
                negChargeTotal += value*(10**pH)/((10**self.aa2chargeNeg[key])+(10**pH)) #using equation provided to return charge from a given residue * its count

        #adding the positive and negative charge of the n and c terminus respectively
        posChargeTotal += (10**self.aaNterm)/((10**self.aaNterm)+(10**pH))
        negChargeTotal += (10**pH)/((10**self.aaCterm)+(10**pH))
        return (posChargeTotal - negChargeTotal)


    def molarExtinction (self):
        """
        Return protein molar extinction coeffecient by summing the molar extinction coeffecients of residues in the protein that absorb at 280.

        input: VLSPADKTNVKAAW
        output: 5500
        """
        molarExtinction = 0 #pre instantiating return value

        for key, value in self.protein.items(): #iterating through protein count dictionary
            if key in ["W","C","Y"]: #searching dictionary of residues that absorb at 280
                molarExtinction += value*self.aa2abs280[key] #multiplying residue count by their extinction ceffecient and summing to the final molar constant
        return molarExtinction


    def massExtinction (self):
        """
        Return mass extinction constant by mutliplying molar extinction coefffecient by molecular weight.

        input: VLSPADKTNVKAAW
        output: 3.67
        """
        myMW =  self.molecularWeight()
        massExtinction = self.molarExtinction() / myMW if myMW else 0.0

        return round(massExtinction,2)


    def molecularWeight (self):
        """
        Return molecular weight by multiplying count of each amino acid by its molecular weight.

        input: VLSPADKTNVKAAW
        output: 1499.7
        """
        sumAA_Weight = 0 #local variable pre instantiation

        for aa in self.washedProteinString: #iterating through washed input containing
            sumAA_Weight += (self.aa2mw[aa]-self.mwH2O)
        sumAA_Weight += self.mwH2O

        return round(sumAA_Weight,1)


