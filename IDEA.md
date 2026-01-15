ORF Detection with Frameshift Identification and Analysis 

Extend ORF analysis by incorporating frameshift detection.

Core functionality:
- Identify the longest ORF in each sequence
- Detect whether the ORF includes a potential **frameshift event**
  
Report:
- The longest ORF containing a frameshift, if present, or
- The longest ORF without a frameshift, if no frameshift is detected.
- Dominance of the longest ORF compared to total ORF in sequence  

Outcome:
- The project should explain the logic used to detect frameshifts and discuss limitations of the approach
- Find the coverage of the longest ORF with a frameshift using bam
- Using sequences in a family do a similarity score between the longest ORF's with a frameshift
- Compare longest ORF to a  database to find closely relative proteins
