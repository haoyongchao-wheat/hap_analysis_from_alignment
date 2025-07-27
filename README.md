# Universal Protein Haplotype Analyzer

A comprehensive tool for analyzing protein haplotypes from amino acid alignment files. 

## Features

- **Multi-format Support**: Accepts FASTA alignment files with automatic sequence processing
- **Flexible Domain Analysis**: Optional domain visualization with customizable colors
- **Publication-Quality Output**: Generates PNG images and editable SVG files optimized for Adobe Illustrator
- **Comprehensive Reports**: Detailed analysis reports with polymorphism and haplotype statistics
- **Cross-Platform**: Available in both Python and R versions
- **Customizable Reference**: Choose any sequence as reference or use the first sequence by default

## Installation

### Python Version

```bash
# Install required dependencies
pip install matplotlib numpy biopython svgwrite
```

### R Version

```r
# Install required packages
install.packages(c("ggplot2", "seqinr", "jsonlite", "svglite", "RColorBrewer"))
```

## Usage

### Python Version

```bash
# Basic usage
python universal_haplotype_analyzer.py input_sequences.fasta

# With domain information
python universal_haplotype_analyzer.py input_sequences.fasta -d domains.txt

# Specify reference sequence and output directory
python universal_haplotype_analyzer.py input_sequences.fasta -r REF_SEQ_ID -o output_folder

# Full command with all options
python universal_haplotype_analyzer.py input_sequences.fasta -d domains.json -r REF_SEQ_ID -o results/
```

### R Version

```bash
# Basic usage
Rscript universal_haplotype_analyzer.R input_sequences.fasta

# With domain information
Rscript universal_haplotype_analyzer.R input_sequences.fasta -d domains.txt

# Specify reference sequence and output directory
Rscript universal_haplotype_analyzer.R input_sequences.fasta -r REF_SEQ_ID -o output_folder
```

## Input File Formats

### FASTA Alignment File

- Standard FASTA format with aligned amino acid sequences
- Stop codons (*) are automatically removed
- Gaps (-) are automatically removed
- First sequence is used as reference by default

**Example:**
```
>NFS10_Reference
MKLLFAILSLVLVAPMAAADTGVLDDFYIKYNGKLINKKTEYEKLSKEYGNKIE...
>Sample_001
MKLLFAILSLVLVAPMAAADTGVLDDFYIKYNGKLINKKTEYEKLSKEYGNKIE...
>Sample_002
MKLLFAILSLVLVAPMAAADTGVLDDFYIKYNGKLINKKTEYEKLSKEYGNKIE...
```

### Domain Information Files

Two formats are supported:

#### 1. Tab-separated Text Format (domains.txt)

```
# Domain file format: name<TAB>start<TAB>end<TAB>[color]
# Lines starting with # are comments
# Positions are 1-based amino acid coordinates
# Color is optional (hex format or color names)

Signal_peptide	1	25	#FF6B6B
Functional_domain	26	120	#4ECDC4
```

#### 2. JSON Format (domains.json)

```json
{
  "Signal_peptide": {
    "start": 1,
    "end": 25,
    "color": "#FF6B6B"
  },
  "Functional_domain": {
    "start": 26,
    "end": 120,
    "color": "#4ECDC4"
  }
}
```

## Output Files

The tool generates several output files:

1. **haplotype_analysis_report.txt**: Comprehensive analysis report
2. **haplotype_diagram.png**: High-resolution PNG image (300 DPI)
3. **haplotype_ai_optimized.svg**: Editable SVG file optimized for Adobe Illustrator

### SVG File Features

- **Layered Structure**: Organized in separate groups for easy editing
- **Unique Element IDs**: Each element has a unique identifier
- **No Overlapping Elements**: Clean structure for professional editing
- **CSS Classes**: Standardized styling for consistent appearance
- **Adobe Illustrator Optimized**: Perfect compatibility with AI for publication-quality figures

## Command Line Options

### Python Version

```
usage: universal_haplotype_analyzer.py [-h] [-d DOMAINS] [-r REFERENCE] [-o OUTPUT] [--version] fasta_file

positional arguments:
  fasta_file            Input FASTA alignment file

optional arguments:
  -h, --help            show this help message and exit
  -d DOMAINS, --domains DOMAINS
                        Domain information file (JSON or tab-separated)
  -r REFERENCE, --reference REFERENCE
                        Reference sequence ID (default: first sequence)
  -o OUTPUT, --output OUTPUT
                        Output directory (default: current directory)
  --version             show program's version number and exit
```

### R Version

```
usage: universal_haplotype_analyzer.R [-h] [-d DOMAINS] [-r REFERENCE] [-o OUTPUT] [--version] fasta_file

Same options as Python version
```

## Examples

### Example 1: Basic Analysis

```bash
python universal_haplotype_analyzer.py final.243.haplotype.prot.fas
```

This will:
- Use the first sequence as reference
- Identify all polymorphic sites
- Classify sequences into haplotypes
- Generate visualization and report

### Example 2: With Domain Information

```bash
python universal_haplotype_analyzer.py final.243.haplotype.prot.fas -d demo_domains.txt
```

This will additionally:
- Display domain boundaries on the diagram
- Color-code domains according to the domain file
- Include domain information in the analysis report

### Example 3: Custom Reference and Output

```bash
python universal_haplotype_analyzer.py final.243.haplotype.prot.fas -r NFS10_Reference -o results/
```

This will:
- Use "NFS10_Reference" as the reference sequence
- Save all output files in the "results/" directory

## Key Features for Publication

### High-Quality Visualizations

- **300 DPI PNG**: Publication-ready raster images
- **Vector SVG**: Scalable graphics for any publication size
- **Professional Color Schemes**: Automatically generated distinct colors
- **Clear Typography**: Arial font for maximum readability

### Adobe Illustrator Optimization

The SVG output is specifically optimized for Adobe Illustrator:

- **Separated Layers**: Background, domains, sequences, polymorphisms, labels, axes, and legend
- **Unique IDs**: Every element can be individually selected and modified
- **Clean Structure**: No overlapping elements that cause editing difficulties
- **Standard CSS**: Consistent styling that can be easily modified

### Comprehensive Analysis

- **Polymorphism Detection**: Identifies all variable sites relative to reference
- **Haplotype Classification**: Groups identical sequences and assigns representative names
- **Statistical Summary**: Polymorphism density and haplotype diversity metrics
- **Detailed Reports**: Complete analysis documentation for methods sections

## Troubleshooting

### Common Issues

1. **"No sequences found"**: Check FASTA file format and ensure sequences are present
2. **"Reference sequence not found"**: Verify the reference ID matches exactly (case-sensitive)
3. **"Domain file error"**: Check domain file format and coordinate ranges
4. **"Permission denied"**: Ensure write permissions for the output directory

### File Format Requirements

- FASTA files must contain aligned amino acid sequences
- Domain coordinates must be 1-based and within sequence length
- JSON domain files must follow the exact structure shown above
- Tab-separated domain files must use actual tab characters (not spaces)

## Citation

If you use this tool in your research, please cite:

```
Universal Protein Haplotype Analyzer v1.0.0
A comprehensive tool for protein haplotype analysis and visualization
```

## License

MIT License - see LICENSE file for details

## Support

For questions, bug reports, or feature requests, please open an issue on GitHub.

## Version History

- **v1.0.0**: Initial release with Python and R versions
  - Multi-format domain support
  - Adobe Illustrator optimized SVG output
  - Comprehensive analysis reports
  - Cross-platform compatibility
