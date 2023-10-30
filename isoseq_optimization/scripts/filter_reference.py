import sys

def filter_fasta_by_length(fasta_file, min_length):
    with open(fasta_file, 'r') as file:
        header = ''
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if len(sequence) >= min_length:
                    print(header)
                    print(sequence)
                header = line
                sequence = ''
            else:
                sequence += line
        if len(sequence) >= min_length:
            print(header)
            print(sequence)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <fasta_file> <minimum_length>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    min_length = int(sys.argv[2])

    filter_fasta_by_length(fasta_file, min_length)
