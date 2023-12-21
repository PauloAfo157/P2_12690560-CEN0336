#!/usr/bin/env python3

import sys
import re

def achar_maior_orf(seq):
    inicio_codon = "ATG"
    fim_codons = ["TAA", "TAG", "TGA"]
    maior_orf = ""

    for frame in range(3):
        for i in range(frame, len(seq), 3):
            codon = seq[i:i+3]
            if codon == inicio_codon:
                orf = ""
                for j in range(i, len(seq), 3):
                    codon = seq[j:j+3]
                    orf += codon
                    if codon in fim_codons:
                        break

                if len(orf) > len(maior_orf):
                    maior_orf = orf

    return maior_orf

def tradu_seq(seq):
    # Dicionário de tradução de códons para aminoácidos
    tabela_codons = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    }

    proteina = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        aminoaci = tabela_codons.get(codon, '')
        if aminoaci == '*':
            break
        proteina += aminoaci

    return proteina

def write_fasta(outfile, identif, seq):
    outfile.write(f">{identif}\n{seq}\n")

def main():
    # Verifica se o número de argumentos é correto
    if len(sys.argv) != 2:
        print("Número de argumentos incorreto")
        sys.exit(1)

    # Obtém o nome do arquivo de entrada a partir da linha de comando
    input_file = sys.argv[1]

    # Define os nomes dos arquivos de saída
    output_fna = "ORF.fna"
    output_faa = "ORF.faa"

    # Abre os arquivos de saída para escrita
    with open(input_file, "r") as infile, open(output_fna, "w") as fna_outfile, open(output_faa, "w") as faa_outfile:
        atual_identif = ""
        atual_seq = ""

        # Loop sobre as linhas do arquivo de entrada
        for l in infile:
            # Verifica o início da sequência
            if l.startswith(">"):
                # Se já existe uma sequência atual
                if atual_identif and atual_seq:
                    process_record(fna_outfile, faa_outfile, atual_identif, atual_seq)

                # Atualiza o identificador e reinicia a sequência
                atual_identif = l.strip()[1:]
                atual_seq = ""
            else:
                # Adiciona a linha na sequência atual
                atual_seq += l.strip()

        # Processa o último registro no arquivo
        if atual_identif and atual_seq:
            process_record(fna_outfile, faa_outfile, atual_identif, atual_seq)

def process_record(fna_outfile, faa_outfile, identif, seq):
    # Encontra o ORF mais longo na sequência
    maior_orf = achar_maior_orf(seq)

    # Traduz o ORF em um peptídeo
    traduzida_seq = tradu_seq(maior_orf)

    # Calcula as coordenadas do ORF
    frame = seq.find(maior_orf) % 3 + 1
    start = seq.find(maior_orf)
    end = start + len(maior_orf)

    # Cria identificadores para os arquivos de saída
    identificador = f"_frame{frame}_{start}_{end}"
    fna_identificador = f"{identif}{identificador}"
    faa_identificador = f"{identif}{identificador}"

    # Escreve as sequências nos arquivos de saída
    write_fasta(fna_outfile, fna_identificador, maior_orf)
    write_fasta(faa_outfile, faa_identificador, traduzida_seq)
