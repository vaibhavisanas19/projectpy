import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO
from Bio import AlignIO, SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

st.title("🧬 Phylogenetic Tree & Conservation Analysis")

st.subheader("📌 Enter FASTA Sequences")
fasta_data = st.text_area(
    "Paste sequences in FASTA format:",
    """>Seq1
ATGCGTACGTTAGTAACTG
>Seq2
ATGCGTACGTTAGTACCTG
>Seq3
ATGCGTACGTTGGTAACTG
>Seq4
ATGCGGACGTTAGTAACTG""",
    height=200,
)


# 🎯 Function to Parse FASTA
def parse_fasta(fasta_text):
    fasta_io = StringIO(fasta_text.strip())
    records = list(SeqIO.parse(fasta_io, "fasta"))
    return records if len(records) >= 2 else None


# 🎯 Function to Build Phylogenetic Tree
def build_phylogenetic_tree(fasta_text):
    fasta_io = StringIO(fasta_text.strip())
    alignment = AlignIO.read(fasta_io, "fasta")
    calculator = DistanceCalculator("identity")
    constructor = DistanceTreeConstructor(calculator, "upgma")
    return constructor.build_tree(alignment)


# 🎯 Function for Conservation Analysis
def conservation_analysis(fasta_text):
    fasta_io = StringIO(fasta_text.strip())
    alignment = AlignIO.read(fasta_io, "fasta")
    conservation_scores_list = []  # Renamed variable to avoid shadowing

    for i in range(alignment.get_alignment_length()):
        column = [seq[i] for seq in alignment]
        most_common = max(set(column), key=column.count)
        score = column.count(most_common) / len(column)
        conservation_scores_list.append(score)  # Using renamed variable

    return conservation_scores_list  # Returning renamed variable


# 🎯 Button to Run Analysis
if st.button("Analyze Sequences"):
    parsed_sequences = parse_fasta(fasta_data)

    if not parsed_sequences:
        st.error("❌ Please provide at least two valid FASTA sequences!")
    else:
        st.success("✅ Phylogenetic Tree Constructed!")
        phylo_tree = build_phylogenetic_tree(fasta_data)

        # Display Tree as ASCII
        st.subheader("Phylogenetic Tree (ASCII)")
        tree_ascii = StringIO()
        Phylo.draw_ascii(phylo_tree, file=tree_ascii)
        st.text(tree_ascii.getvalue())

        # Conservation Analysis
        st.subheader("📊 Conservation Analysis")
        conservation_scores_result = conservation_analysis(fasta_data)  # Using a different variable name

        # Display Conservation Scores
        df_conservation = pd.DataFrame({"Position": range(1, len(conservation_scores_result) + 1),
                                        "Conservation Score": conservation_scores_result})
        st.dataframe(df_conservation)

        # Plot Conservation Scores
        st.subheader("📈 Conservation Score Plot")
        plt.figure(figsize=(8, 4))
        plt.plot(range(1, len(conservation_scores_result) + 1), conservation_scores_result, marker="o", linestyle="-")
        plt.xlabel("Position")
        plt.ylabel("Conservation Score")
        plt.title("Conservation Analysis")
        plt.grid()
        st.pyplot(plt)
