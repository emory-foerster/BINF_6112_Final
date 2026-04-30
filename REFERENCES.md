# References & Acknowledgments

## AI Assistance

This project used Claude (Anthropic) as a development aid during the Spring 2026 semester.

**Citation:**
> Anthropic. (2026). *Claude* (Claude Sonnet 4.6) [Large language model]. https://claude.ai

### What Claude was used for

**General**
- Explaining design decisions and justifying architectural choices
- Advising on modularization strategy and separation of concerns
- Advising on html output formatting

### What Claude was NOT used for
All core biology logic was written by the project team, including:
- ORF detection algorithm (`detect_ORF`, `parse_codons`)
- Frameshift detection logic (`FrameshiftDetector`)
- FASTA file parsing (`fasta_io.py`)
- Gene annotation logic (`get_gene_name`)
- Terminal visualization (`visualize_orf`, `display_gene_coverage`)
- The frameshift sequence track plot (`WIP_report_2.py` core logic)

---

## Academic References

See the Glossary tab of the HTML report (`report2.html`) for the full list of
peer-reviewed references cited in the biological interpretation sections.

Key references:

[1] National Human Genome Research Institute. Open Reading Frame. https://www.genome.gov/genetics-glossary/Open-Reading-Frame

[2] Wu F et al. (2020). A new coronavirus associated with human respiratory disease in China. *Nature* 579:265–269. https://doi.org/10.1038/s41586-020-2008-3

[3] Finkel Y et al. (2021). The coding capacity of SARS-CoV-2. *Nature* 589:125–130. https://doi.org/10.1038/s41586-020-2739-1

[4] Brierley I, Digard P & Inglis SC (1989). Characterization of an efficient coronavirus ribosomal frameshifting signal. *Cell* 57(4):537–547. https://doi.org/10.1016/0092-8674(89)90124-4

[5] Kelly JA et al. (2020). Structural and functional conservation of the programmed -1 ribosomal frameshift signal of SARS-CoV-2. *RNA* 26(9):1175–1189. https://doi.org/10.1261/rna.076141.120

[6] Bhatt PR et al. (2021). Structural basis of ribosomal frameshifting during translation of the SARS-CoV-2 RNA genome. *Science* 372(6548):1306–1313. https://doi.org/10.1126/science.abf3546

[7] Walls AC et al. (2020). Structure, function, and antigenicity of the SARS-CoV-2 spike glycoprotein. *Cell* 181(2):281–292. https://doi.org/10.1016/j.cell.2020.02.058

[8] Cock PJA, et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.

[9] Farabaugh PJ. (1996). Programmed translational frameshifting. Microbiological Reviews, 60(1), 103–134.
