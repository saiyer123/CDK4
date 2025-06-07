from Bio import Entrez, SeqIO, AlignIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
import io
import tempfile


Entrez.email = "l1joseph@ucsd.edu"

def get_complete_protein_sequence(gene_id, species):
    try:
        # First try to get specific curated protein sequences for known orthologs
        protein_ids = {
            "Homo sapiens": "NP_000066.1",        # Human CDK4
            "Pan troglodytes": "XP_016798531.1",  # Chimp CDK4
            "Mus musculus": "NP_034007.3",        # Mouse Cdk4
            "Danio rerio": "NP_571420.1",         # Zebrafish cdk4
            "Drosophila melanogaster": "NP_523355.2", # Fruit fly Cdk4
            "Caenorhabditis elegans": "NP_001022846.1", # Worm cdk-4
            "Saccharomyces cerevisiae": "NP_009718.3"  # Yeast CDC28
        }
        
        if species in protein_ids:
            try:
                handle = Entrez.efetch(db="protein", id=protein_ids[species], rettype="fasta", retmode="text")
                sequence_record = SeqIO.read(handle, "fasta")
                handle.close()
                print(f"Retrieved curated {species} sequence: {len(sequence_record.seq)} aa")
                return sequence_record
            except Exception as e:
                print(f"Could not retrieve curated sequence for {species}: {e}")
                # Fall back to gene-based lookup
        
        # Get protein IDs linked to this gene
        handle = Entrez.elink(dbfrom="gene", db="protein", id=gene_id)
        record = Entrez.read(handle)
        handle.close()
        
        if record[0]["LinkSetDb"] and len(record[0]["LinkSetDb"][0]["Link"]) > 0:
            # Get all protein IDs
            protein_ids = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
            
            # Try to find the longest protein sequence (likely the most complete)
            longest_seq = None
            max_length = 0
            
            for protein_id in protein_ids[:3]:  # Check first 3 to avoid too many queries
                handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(handle, "fasta")
                handle.close()
                
                if len(seq_record.seq) > max_length:
                    max_length = len(seq_record.seq)
                    longest_seq = seq_record
                
                time.sleep(1)  # Be nice to NCBI servers
            
            if longest_seq:
                return longest_seq
            else:
                print(f"No suitable protein sequence found for {species}")
                return None
        else:
            print(f"No protein links found for gene {gene_id} in {species}")
            return None
    
    except Exception as e:
        print(f"Error retrieving protein sequence for {species}: {e}")
        return None

def perform_multiple_sequence_alignment(sequences):
    # Create a temporary file for the sequences
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fasta") as temp_file:
        for species, seq_record in sequences.items():
            # Use shorter names for better display
            short_name = f"{species.split()[0]}_{species.split()[1][:3]}" if " " in species else species
            seq_record.id = short_name
            seq_record.description = ""
            SeqIO.write(seq_record, temp_file, "fasta")
        
        temp_filename = temp_file.name
    
    # Use human CDK4 as the reference sequence
    reference_seq = sequences["Homo sapiens"]
    
    # Perform BLASTp alignment
    try:
        print("Running BLASTp alignment...")
        result_handle = NCBIWWW.qblast("blastp", "nr", reference_seq.seq)
        blast_records = NCBIXML.parse(result_handle)
        
        # Create alignment from BLAST results
        aligned_sequences = []
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Create a sequence record for the aligned sequence
                    aligned_seq = SeqIO.SeqRecord(
                        seq=hsp.sbjct,
                        id=alignment.title.split()[0],
                        description=""
                    )
                    aligned_sequences.append(aligned_seq)
        
        # Create a multiple sequence alignment
        alignment = AlignIO.MultipleSeqAlignment(aligned_sequences)
        print(f"BLASTp alignment completed with {len(alignment)} sequences")
        
        return alignment
    except Exception as e:
        print(f"Error performing BLASTp alignment: {e}")
        print("Falling back to pairwise alignments")
        return None

def analyze_cdk4_evolution_with_phylogeny():    
    # List of key species to search for orthologs
    key_species = [
        "Homo sapiens", "Pan troglodytes", "Mus musculus", 
        "Gallus gallus", "Danio rerio", "Drosophila melanogaster",
        "Caenorhabditis elegans", "Saccharomyces cerevisiae"
    ]
    
    # Estimated divergence times (millions of years ago) for reference
    divergence_times = {
        "Homo sapiens": 0,
        "Pan troglodytes": 6,
        "Mus musculus": 90,
        "Gallus gallus": 310,
        "Danio rerio": 435,
        "Drosophila melanogaster": 780,
        "Caenorhabditis elegans": 900,
        "Saccharomyces cerevisiae": 1500
    }
    
    # Get ortholog sequences
    sequences = {}
    orthologs = []
    
    for species in key_species:
        try:
            # Search for CDK4 in the species
            search_term = f"CDK4[Gene Name] AND {species}[Organism]"
            if species == "Saccharomyces cerevisiae":
                # Yeast CDK4 ortholog is actually CDC28
                search_term = f"CDC28[Gene Name] AND {species}[Organism]"
            elif species == "Drosophila melanogaster":
                # Fruit fly ortholog is Cdk4/Cdk4-like
                search_term = f"Cdk4[Gene Name] AND {species}[Organism]"
            elif species == "Caenorhabditis elegans":
                # C. elegans ortholog is cdk-4
                search_term = f"cdk-4[Gene Name] AND {species}[Organism]"
                
            handle = Entrez.esearch(db="gene", term=search_term)
            record = Entrez.read(handle)
            handle.close()
            
            if len(record["IdList"]) > 0:
                gene_id = record["IdList"][0]
                
                # Get gene details
                handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
                gene_record = Entrez.read(handle)
                handle.close()
                
                if len(gene_record) > 0:
                    gene_data = gene_record[0]
                    gene_symbol = gene_data.get("Entrezgene_gene", {}).get("Gene-ref", {}).get("Gene-ref_locus", "")
                    
                    if not gene_symbol and "Entrezgene_gene" in gene_data:
                        gene_symbol = gene_data["Entrezgene_gene"]["Gene-ref"].get("Gene-ref_locus_tag", "")
                    
                    if not gene_symbol:
                        gene_symbol = search_term.split("[")[0]
                    
                    # Get complete protein sequence
                    print(f"Retrieving complete protein sequence for {gene_symbol} in {species}...")
                    sequence = get_complete_protein_sequence(gene_id, species)
                    
                    if sequence:
                        sequences[species] = sequence
                        orthologs.append({
                            "Species": species,
                            "Gene_Symbol": gene_symbol,
                            "NCBI_Gene_ID": gene_id,
                            "Sequence_Length": len(sequence.seq),
                            "Known_Divergence_MYA": divergence_times.get(species, "Unknown")
                        })
                        
                        print(f"Found {gene_symbol} in {species}, sequence length: {len(sequence.seq)}")
                    
                    # # Avoid overwhelming NCBI servers
                    # time.sleep(1)
            
        except Exception as e:
            print(f"Error retrieving {species} ortholog: {e}")
    
    # Create DataFrame
    orthologs_df = pd.DataFrame(orthologs)
    
    if not orthologs_df.empty and len(sequences) >= 3:
        print("\nCDK4 Orthologs:")
        print(orthologs_df[["Species", "Gene_Symbol", "Sequence_Length", "Known_Divergence_MYA"]].to_string(index=False))
        
        # Perform multiple sequence alignment
        alignment = perform_multiple_sequence_alignment(sequences)
        
        if alignment:
            # Calculate distance matrix
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(alignment)
            
            # Build phylogenetic tree
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)
            
            # Display the tree
            plt.figure(figsize=(10, 8))
            Phylo.draw(tree, do_show=False)
            plt.title("Phylogenetic Tree of CDK4 Orthologs")
            plt.savefig("cdk4_phylogenetic_tree.png", dpi=300, bbox_inches='tight')
            plt.show()
            
            # Find most distant ortholog
            species_names = [seq.id for seq in alignment]
            human_index = next((i for i, name in enumerate(species_names) if name.startswith("Homo")), 0)
            
            distances = []
            for i, name in enumerate(species_names):
                if i != human_index:
                    species = next((s for s in key_species if name.startswith(s.split()[0])), name)
                    gene_symbol = orthologs_df[orthologs_df["Species"].str.startswith(species.split()[0])]["Gene_Symbol"].values[0]
                    distances.append((species, gene_symbol, dm[human_index, i]))
            
            # Sort by distance
            distances.sort(key=lambda x: x[2], reverse=True)
            
            print("\nEvolutionary Distances from Human CDK4 (based on multiple alignment):")
            for species, gene, distance in distances:
                print(f"{species} ({gene}): {distance:.4f}")
            
            most_distant = distances[0]
            print(f"\nMost distant ortholog from human CDK4: {most_distant[1]} in {most_distant[0]} "
                  f"(distance: {most_distant[2]:.4f})")
            
            # Get known divergence time for most distant species
            most_distant_species = most_distant[0]
            most_distant_time = next((divergence_times.get(s, "Unknown") for s in key_species if s.startswith(most_distant_species.split()[0])), "Unknown")
            
            
            if most_distant_time != "Unknown":
                print(f"{most_distant_time}, {most_distant[1]}, {most_distant_species}")
            
            return orthologs_df, sequences, alignment, tree
        
    return orthologs_df, sequences, None, None


orthologs_df, sequences, alignment, tree = analyze_cdk4_evolution_with_phylogeny()
