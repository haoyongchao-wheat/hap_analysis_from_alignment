#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Universal Protein Haplotype Analyzer

A comprehensive tool for analyzing protein haplotypes from amino acid alignment files.
Generates high-quality visualizations and detailed reports suitable for publication.

Author: Haplotype Analysis Tool
Version: 1.0.0
License: MIT
"""

import argparse
import os
import sys
from pathlib import Path
import json
from collections import defaultdict, Counter
import random
import colorsys

# Core libraries
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from Bio import SeqIO
import svgwrite

# Set matplotlib parameters
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.unicode_minus'] = False

class UniversalHaplotypeAnalyzer:
    """
    Universal Protein Haplotype Analyzer
    
    Analyzes protein haplotypes from amino acid alignment files and generates
    publication-quality visualizations and comprehensive reports.
    """
    
    def __init__(self, fasta_file, domain_file=None, reference_id=None, output_dir=None):
        """
        Initialize the analyzer
        
        Args:
            fasta_file (str): Path to the FASTA alignment file
            domain_file (str, optional): Path to the domain information file
            reference_id (str, optional): ID of reference sequence (default: first sequence)
            output_dir (str, optional): Output directory (default: current directory)
        """
        self.fasta_file = Path(fasta_file)
        self.domain_file = Path(domain_file) if domain_file else None
        self.reference_id = reference_id
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()
        
        # Data storage
        self.sequences = {}
        self.reference_seq = None
        self.reference_name = None
        self.haplotypes = {}
        self.polymorphisms = []
        self.domains = {}
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate input files
        self._validate_inputs()
    
    def _validate_inputs(self):
        """Validate input files"""
        if not self.fasta_file.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_file}")
        
        if self.domain_file and not self.domain_file.exists():
            raise FileNotFoundError(f"Domain file not found: {self.domain_file}")
    
    def load_sequences(self):
        """Load sequences from FASTA file"""
        print(f"Loading sequences from {self.fasta_file}...")
        
        sequence_order = []
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            # Remove stop codons (*) and gaps (-) but keep other characters
            seq_str = str(record.seq).replace('*', '').replace('-', '')
            self.sequences[record.id] = seq_str
            sequence_order.append(record.id)
            
        if not self.sequences:
            raise ValueError("No sequences found in FASTA file")
        
        # Set reference sequence
        if self.reference_id:
            if self.reference_id not in self.sequences:
                raise ValueError(f"Reference sequence '{self.reference_id}' not found")
            self.reference_name = self.reference_id
        else:
            # Use first sequence as reference
            self.reference_name = sequence_order[0]
        
        self.reference_seq = self.sequences[self.reference_name]
        
        print(f"Loaded {len(self.sequences)} sequences")
        print(f"Reference sequence: {self.reference_name} (length: {len(self.reference_seq)} aa)")
    
    def load_domains(self):
        """Load domain information from file"""
        if not self.domain_file:
            print("No domain file provided, skipping domain analysis")
            return
        
        print(f"Loading domain information from {self.domain_file}...")
        
        try:
            if self.domain_file.suffix.lower() == '.json':
                self._load_domains_json()
            else:
                self._load_domains_text()
        except Exception as e:
            print(f"Warning: Could not load domain file: {e}")
            print("Continuing without domain information")
            self.domains = {}
    
    def _load_domains_json(self):
        """Load domains from JSON file"""
        with open(self.domain_file, 'r') as f:
            domain_data = json.load(f)
        
        for domain_name, info in domain_data.items():
            start = info['start']
            end = info['end']
            color = info.get('color', None)
            self.domains[domain_name] = (start, end, color)
    
    def _load_domains_text(self):
        """Load domains from text file (tab-separated: name start end [color])"""
        with open(self.domain_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split('\t')
                if len(parts) < 3:
                    print(f"Warning: Invalid line {line_num} in domain file: {line}")
                    continue
                
                domain_name = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    color = parts[3] if len(parts) > 3 else None
                    self.domains[domain_name] = (start, end, color)
                except ValueError:
                    print(f"Warning: Invalid coordinates in line {line_num}: {line}")
    
    def _generate_colors(self, n_colors):
        """Generate n distinct colors using HSV color space"""
        colors = []
        for i in range(n_colors):
            hue = i / n_colors
            saturation = 0.7 + (i % 2) * 0.2  # Alternate between 0.7 and 0.9
            value = 0.8 + (i % 3) * 0.1       # Vary brightness slightly
            rgb = colorsys.hsv_to_rgb(hue, saturation, value)
            hex_color = '#{:02x}{:02x}{:02x}'.format(
                int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255)
            )
            colors.append(hex_color)
        return colors
    
    def _assign_domain_colors(self):
        """Assign colors to domains"""
        if not self.domains:
            return
        
        # Check if any domains already have colors
        domains_without_colors = [name for name, (start, end, color) in self.domains.items() if not color]
        
        if domains_without_colors:
            # Generate colors for domains without colors
            new_colors = self._generate_colors(len(domains_without_colors))
            
            for i, domain_name in enumerate(domains_without_colors):
                start, end, _ = self.domains[domain_name]
                self.domains[domain_name] = (start, end, new_colors[i])
        
        print(f"Assigned colors to {len(self.domains)} domains")
    
    def find_polymorphisms(self):
        """Find polymorphic sites relative to reference sequence"""
        print("Identifying polymorphic sites...")
        
        ref_len = len(self.reference_seq)
        
        for pos in range(ref_len):
            ref_aa = self.reference_seq[pos]
            variants = set()
            
            for seq_id, sequence in self.sequences.items():
                if seq_id == self.reference_name:
                    continue
                
                if pos < len(sequence):
                    if sequence[pos] != ref_aa:
                        variants.add(sequence[pos])
            
            if variants:
                self.polymorphisms.append({
                    'position': pos + 1,  # 1-based position
                    'reference': ref_aa,
                    'variants': list(variants)
                })
        
        print(f"Found {len(self.polymorphisms)} polymorphic sites")
    
    def classify_haplotypes(self):
        """Classify and merge haplotypes"""
        print("Classifying haplotypes...")
        
        haplotype_groups = defaultdict(list)
        
        for seq_id, sequence in self.sequences.items():
            if seq_id == self.reference_name:
                continue
            
            # Generate haplotype signature based on polymorphic sites
            haplotype_signature = []
            for poly in self.polymorphisms:
                pos = poly['position'] - 1  # Convert to 0-based
                if pos < len(sequence):
                    haplotype_signature.append(f"{poly['position']}:{sequence[pos]}")
                else:
                    haplotype_signature.append(f"{poly['position']}:-")
            
            signature_str = '|'.join(haplotype_signature)
            haplotype_groups[signature_str].append(seq_id)
        
        # Assign representative names to haplotype groups
        haplotype_counter = 1
        for signature, seq_ids in haplotype_groups.items():
            if len(seq_ids) == 1:
                # Single sequence - use its name
                haplotype_name = seq_ids[0]
            else:
                # Multiple sequences - create group name
                haplotype_name = f"Haplotype_{haplotype_counter}"
                haplotype_counter += 1
            
            # Use first sequence as representative
            representative_seq = self.sequences[seq_ids[0]]
            self.haplotypes[haplotype_name] = {
                'sequence': representative_seq,
                'members': seq_ids,
                'signature': signature
            }
        
        print(f"Identified {len(self.haplotypes)} distinct haplotypes")
        for name, info in self.haplotypes.items():
            print(f"  {name}: {len(info['members'])} member(s)")
    
    def generate_report(self):
        """Generate comprehensive analysis report"""
        report_file = self.output_dir / "haplotype_analysis_report.txt"
        
        print(f"Generating analysis report: {report_file}")
        
        with open(report_file, 'w') as f:
            f.write("PROTEIN HAPLOTYPE ANALYSIS REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            # Basic information
            f.write("BASIC INFORMATION\n")
            f.write("-" * 20 + "\n")
            f.write(f"Input file: {self.fasta_file.name}\n")
            f.write(f"Reference sequence: {self.reference_name}\n")
            f.write(f"Reference length: {len(self.reference_seq)} amino acids\n")
            f.write(f"Total sequences analyzed: {len(self.sequences)}\n")
            f.write(f"Domain file: {self.domain_file.name if self.domain_file else 'None'}\n")
            f.write(f"Number of domains: {len(self.domains)}\n\n")
            
            # Polymorphism summary
            f.write("POLYMORPHISM ANALYSIS\n")
            f.write("-" * 25 + "\n")
            f.write(f"Total polymorphic sites: {len(self.polymorphisms)}\n")
            f.write(f"Polymorphism density: {len(self.polymorphisms)/len(self.reference_seq)*100:.2f}%\n\n")
            
            if self.polymorphisms:
                f.write("Polymorphic sites details:\n")
                for i, poly in enumerate(self.polymorphisms, 1):
                    f.write(f"  {i:2d}. Position {poly['position']:4d}: {poly['reference']} -> {', '.join(poly['variants'])}\n")
                f.write("\n")
            
            # Haplotype summary
            f.write("HAPLOTYPE ANALYSIS\n")
            f.write("-" * 20 + "\n")
            f.write(f"Number of distinct haplotypes: {len(self.haplotypes)}\n\n")
            
            for i, (name, info) in enumerate(self.haplotypes.items(), 1):
                f.write(f"Haplotype {i}: {name}\n")
                f.write(f"  Sequence length: {len(info['sequence'])} amino acids\n")
                f.write(f"  Number of members: {len(info['members'])}\n")
                f.write(f"  Members: {', '.join(info['members'])}\n")
                
                # Calculate differences from reference
                differences = 0
                for poly in self.polymorphisms:
                    pos = poly['position'] - 1
                    if pos < len(info['sequence']) and pos < len(self.reference_seq):
                        if info['sequence'][pos] != self.reference_seq[pos]:
                            differences += 1
                
                f.write(f"  Differences from reference: {differences} sites\n")
                
                # Check for truncation
                if len(info['sequence']) < len(self.reference_seq):
                    f.write(f"  Status: Truncated at position {len(info['sequence'])}\n")
                else:
                    f.write(f"  Status: Full length\n")
                f.write("\n")
            
            # Domain information
            if self.domains:
                f.write("DOMAIN INFORMATION\n")
                f.write("-" * 18 + "\n")
                for domain_name, (start, end, color) in self.domains.items():
                    f.write(f"  {domain_name}: {start}-{end} aa (color: {color})\n")
                f.write("\n")
            
            # Sequence statistics
            f.write("SEQUENCE STATISTICS\n")
            f.write("-" * 19 + "\n")
            lengths = [len(seq) for seq in self.sequences.values()]
            f.write(f"  Mean length: {np.mean(lengths):.1f} aa\n")
            f.write(f"  Length range: {min(lengths)}-{max(lengths)} aa\n")
            f.write(f"  Standard deviation: {np.std(lengths):.1f} aa\n")
            
            # Length distribution
            length_counts = Counter(lengths)
            f.write(f"\n  Length distribution:\n")
            for length in sorted(length_counts.keys()):
                count = length_counts[length]
                f.write(f"    {length} aa: {count} sequence(s)\n")
        
        return report_file
    
    def create_png_visualization(self):
        """Create PNG visualization"""
        png_file = self.output_dir / "haplotype_diagram.png"
        
        print(f"Creating PNG visualization: {png_file}")
        
        # Calculate figure dimensions
        num_sequences = len(self.haplotypes) + 1
        seq_length = len(self.reference_seq)
        
        fig_width = 16
        fig_height = max(8, num_sequences * 0.8)
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)
        
        # Draw sequences
        y_pos = num_sequences - 1
        
        # Reference sequence
        self._draw_sequence_png(ax, y_pos, self.reference_name, self.reference_seq, 
                               seq_length, is_reference=True)
        y_pos -= 1
        
        # Haplotype sequences
        for haplotype_name, info in self.haplotypes.items():
            display_name = haplotype_name
            if len(info['members']) > 1:
                display_name += f" ({len(info['members'])} members)"
            
            self._draw_sequence_png(ax, y_pos, display_name, info['sequence'], 
                                   seq_length, is_reference=False)
            y_pos -= 1
        
        # Setup axes
        self._setup_png_axes(ax, seq_length, num_sequences)
        
        # Add legend
        if self.domains:
            self._add_png_legend(ax, seq_length)
        
        plt.tight_layout()
        plt.savefig(png_file, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        return png_file
    
    def _draw_sequence_png(self, ax, y_pos, name, sequence, max_length, is_reference=False):
        """Draw a single sequence in PNG format"""
        actual_length = len(sequence)
        bar_height = 0.6
        
        # Sequence bar
        if is_reference:
            bar_color = '#ECF0F1'
            edge_color = '#2C3E50'
        else:
            bar_color = '#FFFFFF'
            edge_color = '#2C3E50'
        
        rect = patches.Rectangle((0, y_pos - bar_height/2), actual_length, bar_height,
                               linewidth=2, edgecolor=edge_color, facecolor=bar_color)
        ax.add_patch(rect)
        
        # Domains
        if self.domains:
            for domain_name, (start, end, color) in self.domains.items():
                if start <= actual_length:
                    domain_start = max(0, start - 1)
                    domain_end = min(end, actual_length)
                    if domain_end > domain_start:
                        domain_rect = patches.Rectangle(
                            (domain_start, y_pos - bar_height/2), 
                            domain_end - domain_start, bar_height,
                            linewidth=1, edgecolor='white', 
                            facecolor=color, alpha=0.8
                        )
                        ax.add_patch(domain_rect)
        
        # Truncation markers
        if actual_length < max_length:
            ax.plot([actual_length, max_length], [y_pos, y_pos], 
                   linestyle='--', color='#95A5A6', linewidth=2, alpha=0.7)
            ax.plot(actual_length, y_pos, marker='|', markersize=15, 
                   color='#E74C3C', markeredgewidth=4)
        
        # Polymorphisms
        if not is_reference:
            for poly in self.polymorphisms:
                pos = poly['position'] - 1
                if pos < actual_length and pos < len(self.reference_seq):
                    if sequence[pos] != self.reference_seq[pos]:
                        ax.plot([pos, pos], [y_pos - 0.4, y_pos + 0.4], 
                               color='#E74C3C', linewidth=3, alpha=0.9)
        
        # Label
        ax.text(-max_length*0.05, y_pos, name, ha='right', va='center', 
               fontsize=12, fontweight='bold')
    
    def _setup_png_axes(self, ax, seq_length, num_sequences):
        """Setup axes for PNG visualization"""
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['bottom'].set_linewidth(2)
        
        ax.set_xlim(0, seq_length)
        ax.set_ylim(-0.8, num_sequences - 0.2)
        ax.set_xlabel('Amino acid position', fontsize=14, fontweight='bold')
        ax.set_ylabel('')
        
        x_ticks = np.arange(0, seq_length + 1, max(100, seq_length // 10))
        ax.set_xticks(x_ticks)
        ax.tick_params(axis='x', which='major', labelsize=12, length=8, width=2)
        ax.tick_params(axis='y', which='both', left=False, labelleft=False)
    
    def _add_png_legend(self, ax, seq_length):
        """Add legend to PNG visualization"""
        legend_elements = []
        
        for domain_name, (start, end, color) in self.domains.items():
            legend_elements.append(
                patches.Patch(color=color, alpha=0.8, 
                            label=f'{domain_name} ({start}-{end} aa)')
            )
        
        # Add polymorphism legend
        legend_elements.append(
            plt.Line2D([0], [0], color='#E74C3C', linewidth=3, 
                      label='Polymorphic sites')
        )
        
        ax.legend(handles=legend_elements, loc='upper right', 
                 bbox_to_anchor=(1, 1), fontsize=10)
    
    def create_svg_visualization(self):
        """Create editable SVG visualization"""
        svg_file = self.output_dir / "haplotype_diagram_editable.svg"
        
        print(f"Creating editable SVG visualization: {svg_file}")
        
        # Calculate dimensions
        num_sequences = len(self.haplotypes) + 1
        seq_length = len(self.reference_seq)
        
        # SVG canvas dimensions
        margin = 150
        canvas_width = 1600
        canvas_height = max(800, num_sequences * 80 + 200)
        
        # Drawing area
        plot_width = canvas_width - 2 * margin
        scale_x = plot_width / seq_length
        row_height = 60
        
        # Create SVG document
        dwg = svgwrite.Drawing(str(svg_file), size=(canvas_width, canvas_height))
        
        # Add CSS styles
        self._add_svg_styles(dwg)
        
        # Create layer groups
        main_group = dwg.g(id="haplotype_diagram")
        sequences_bg_group = dwg.g(id="sequences_background")
        domains_group = dwg.g(id="protein_domains")
        sequences_outline_group = dwg.g(id="sequences_outline")
        polymorphisms_group = dwg.g(id="polymorphic_sites")
        truncation_group = dwg.g(id="truncation_markers")
        labels_group = dwg.g(id="sequence_labels")
        axis_group = dwg.g(id="coordinate_axis")
        legend_group = dwg.g(id="figure_legend")
        
        # Draw sequences
        y_start = 100
        current_y = y_start
        
        # Reference sequence
        self._draw_sequence_svg(dwg, sequences_bg_group, domains_group, 
                               sequences_outline_group, polymorphisms_group,
                               truncation_group, labels_group,
                               "reference", self.reference_name, self.reference_seq,
                               current_y, margin, scale_x, seq_length, True)
        current_y += row_height
        
        # Haplotype sequences
        haplotype_index = 0
        for haplotype_name, info in self.haplotypes.items():
            display_name = haplotype_name
            if len(info['members']) > 1:
                display_name += f" ({len(info['members'])} members)"
            
            self._draw_sequence_svg(dwg, sequences_bg_group, domains_group,
                                   sequences_outline_group, polymorphisms_group,
                                   truncation_group, labels_group,
                                   f"haplotype_{haplotype_index}", display_name, 
                                   info['sequence'], current_y, margin, scale_x, 
                                   seq_length, False)
            current_y += row_height
            haplotype_index += 1
        
        # Add coordinate axis
        self._add_svg_axis(dwg, axis_group, current_y, margin, scale_x, seq_length)
        
        # Add legend
        if self.domains:
            self._add_svg_legend(dwg, legend_group, canvas_width)
        
        # Assemble all groups
        main_group.add(sequences_bg_group)
        main_group.add(domains_group)
        main_group.add(sequences_outline_group)
        main_group.add(polymorphisms_group)
        main_group.add(truncation_group)
        main_group.add(labels_group)
        main_group.add(axis_group)
        main_group.add(legend_group)
        
        dwg.add(main_group)
        dwg.save()
        
        return svg_file
    
    def _add_svg_styles(self, dwg):
        """Add CSS styles to SVG"""
        dwg.defs.add(dwg.style("""
            .sequence-outline { fill: none; stroke: #2C3E50; stroke-width: 2; }
            .reference-fill { fill: #ECF0F1; }
            .haplotype-fill { fill: #FFFFFF; }
            .polymorphism-line { stroke: #E74C3C; stroke-width: 3; fill: none; }
            .truncation-line { stroke: #95A5A6; stroke-width: 2; stroke-dasharray: 8,4; fill: none; }
            .termination-marker { stroke: #E74C3C; stroke-width: 4; fill: none; }
            .sequence-label { font-family: Arial, sans-serif; font-size: 14px; font-weight: bold; fill: #2C3E50; }
            .axis-label { font-family: Arial, sans-serif; font-size: 12px; fill: #34495E; }
            .axis-title { font-family: Arial, sans-serif; font-size: 14px; font-weight: bold; fill: #2C3E50; }
            .legend-text { font-family: Arial, sans-serif; font-size: 11px; fill: #2C3E50; }
        """))
    
    def _draw_sequence_svg(self, dwg, bg_group, domains_group, outline_group, 
                          poly_group, trunc_group, labels_group,
                          seq_id, display_name, sequence, y_pos, margin, 
                          scale_x, max_length, is_reference):
        """Draw a single sequence in SVG format"""
        actual_length = len(sequence)
        bar_height = 50
        
        # Background rectangle
        bg_class = "reference-fill" if is_reference else "haplotype-fill"
        bg_rect = dwg.rect(
            insert=(margin, y_pos - bar_height//2),
            size=(actual_length * scale_x, bar_height),
            class_=bg_class,
            id=f"{seq_id}_background"
        )
        bg_group.add(bg_rect)
        
        # Outline
        outline_rect = dwg.rect(
            insert=(margin, y_pos - bar_height//2),
            size=(actual_length * scale_x, bar_height),
            class_="sequence-outline",
            id=f"{seq_id}_outline"
        )
        outline_group.add(outline_rect)
        
        # Domains
        if self.domains:
            for domain_name, (start, end, color) in self.domains.items():
                if start <= actual_length:
                    domain_start = max(0, start - 1)
                    domain_end = min(end, actual_length)
                    if domain_end > domain_start:
                        domain_rect = dwg.rect(
                            insert=(margin + domain_start * scale_x, y_pos - bar_height//2),
                            size=((domain_end - domain_start) * scale_x, bar_height),
                            fill=color,
                            stroke="white",
                            stroke_width=1,
                            opacity=0.8,
                            id=f"{seq_id}_domain_{domain_name.lower().replace(' ', '_')}"
                        )
                        domains_group.add(domain_rect)
        
        # Truncation markers
        if actual_length < max_length:
            # Dashed line
            trunc_line = dwg.line(
                start=(margin + actual_length * scale_x, y_pos),
                end=(margin + max_length * scale_x, y_pos),
                class_="truncation-line",
                id=f"{seq_id}_truncation_line"
            )
            trunc_group.add(trunc_line)
            
            # Termination marker
            term_marker = dwg.line(
                start=(margin + actual_length * scale_x, y_pos - 20),
                end=(margin + actual_length * scale_x, y_pos + 20),
                class_="termination-marker",
                id=f"{seq_id}_termination_marker"
            )
            trunc_group.add(term_marker)
        
        # Polymorphisms
        if not is_reference:
            poly_index = 0
            for poly in self.polymorphisms:
                pos = poly['position'] - 1
                if pos < actual_length and pos < len(self.reference_seq):
                    if sequence[pos] != self.reference_seq[pos]:
                        poly_line = dwg.line(
                            start=(margin + pos * scale_x, y_pos - 30),
                            end=(margin + pos * scale_x, y_pos + 30),
                            class_="polymorphism-line",
                            id=f"{seq_id}_polymorphism_{poly_index}"
                        )
                        poly_group.add(poly_line)
                        poly_index += 1
        
        # Label
        label = dwg.text(display_name,
                        insert=(margin - 10, y_pos + 5),
                        class_="sequence-label",
                        text_anchor="end",
                        id=f"{seq_id}_label")
        labels_group.add(label)
    
    def _add_svg_axis(self, dwg, axis_group, y_pos, margin, scale_x, seq_length):
        """Add coordinate axis to SVG"""
        axis_y = y_pos + 20
        
        # Main axis line
        axis_line = dwg.line(
            start=(margin, axis_y),
            end=(margin + seq_length * scale_x, axis_y),
            stroke="#2C3E50",
            stroke_width=2,
            id="x_axis_line"
        )
        axis_group.add(axis_line)
        
        # Ticks and labels
        tick_interval = max(100, seq_length // 10)
        for tick in range(0, seq_length + 1, tick_interval):
            tick_x = margin + tick * scale_x
            
            # Tick line
            tick_line = dwg.line(
                start=(tick_x, axis_y),
                end=(tick_x, axis_y + 8),
                stroke="#2C3E50",
                stroke_width=2,
                id=f"x_tick_{tick}"
            )
            axis_group.add(tick_line)
            
            # Tick label
            tick_label = dwg.text(str(tick),
                                 insert=(tick_x, axis_y + 25),
                                 class_="axis-label",
                                 text_anchor="middle",
                                 id=f"x_label_{tick}")
            axis_group.add(tick_label)
        
        # Axis title
        axis_title = dwg.text('Amino acid position',
                             insert=(margin + seq_length * scale_x / 2, axis_y + 50),
                             class_="axis-title",
                             text_anchor="middle",
                             id="x_axis_title")
        axis_group.add(axis_title)
    
    def _add_svg_legend(self, dwg, legend_group, canvas_width):
        """Add legend to SVG"""
        legend_x = canvas_width - 350
        legend_y = 100
        
        # Legend title
        legend_title = dwg.text('Protein Domains',
                               insert=(legend_x, legend_y - 10),
                               class_="axis-title",
                               id="legend_title")
        legend_group.add(legend_title)
        
        # Domain legend items
        for i, (domain_name, (start, end, color)) in enumerate(self.domains.items()):
            legend_rect = dwg.rect(
                insert=(legend_x, legend_y + i * 25),
                size=(25, 18),
                fill=color,
                stroke="#2C3E50",
                stroke_width=1,
                id=f"legend_rect_{domain_name.lower().replace(' ', '_')}"
            )
            legend_group.add(legend_rect)
            
            legend_text = dwg.text(f'{domain_name} ({start}-{end} aa)',
                                  insert=(legend_x + 35, legend_y + i * 25 + 14),
                                  class_="legend-text",
                                  id=f"legend_text_{domain_name.lower().replace(' ', '_')}")
            legend_group.add(legend_text)
        
        # Polymorphism legend
        poly_y = legend_y + len(self.domains) * 25 + 15
        poly_line = dwg.line(
            start=(legend_x, poly_y),
            end=(legend_x + 25, poly_y),
            class_="polymorphism-line",
            id="legend_polymorphism_line"
        )
        legend_group.add(poly_line)
        
        poly_text = dwg.text('Polymorphic sites',
                            insert=(legend_x + 35, poly_y + 5),
                            class_="legend-text",
                            id="legend_polymorphism_text")
        legend_group.add(poly_text)
    
    def run_analysis(self):
        """Run complete haplotype analysis"""
        print("Starting Universal Protein Haplotype Analysis...")
        print("=" * 50)
        
        try:
            # Load data
            self.load_sequences()
            self.load_domains()
            self._assign_domain_colors()
            
            # Analyze
            self.find_polymorphisms()
            self.classify_haplotypes()
            
            # Generate outputs
            report_file = self.generate_report()
            png_file = self.create_png_visualization()
            svg_file = self.create_svg_visualization()
            
            print("\n" + "=" * 50)
            print("Analysis completed successfully!")
            print("\nGenerated files:")
            print(f"  - Report: {report_file}")
            print(f"  - PNG visualization: {png_file}")
            print(f"  - Editable SVG: {svg_file}")
            
            return {
                'report': report_file,
                'png': png_file,
                'svg': svg_file
            }
            
        except Exception as e:
            print(f"\nError during analysis: {e}")
            raise

def main():
    """Main function with command line interface"""
    parser = argparse.ArgumentParser(
        description='Universal Protein Haplotype Analyzer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis
  python universal_haplotype_analyzer.py sequences.fasta
  
  # With domain information
  python universal_haplotype_analyzer.py sequences.fasta -d domains.txt
  
  # Specify reference sequence and output directory
  python universal_haplotype_analyzer.py sequences.fasta -r REF_SEQ -o results/
        """
    )
    
    parser.add_argument('fasta_file', 
                       help='Input FASTA file with aligned amino acid sequences')
    parser.add_argument('-d', '--domains', 
                       help='Domain information file (tab-separated or JSON)')
    parser.add_argument('-r', '--reference', 
                       help='Reference sequence ID (default: first sequence)')
    parser.add_argument('-o', '--output', 
                       help='Output directory (default: current directory)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
    
    args = parser.parse_args()
    
    try:
        analyzer = UniversalHaplotypeAnalyzer(
            fasta_file=args.fasta_file,
            domain_file=args.domains,
            reference_id=args.reference,
            output_dir=args.output
        )
        
        analyzer.run_analysis()
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()