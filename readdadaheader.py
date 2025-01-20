import sys
from baseband import dada

def print_dada_header(filepath):
    """Print the header of a DADA file."""
    with dada.open(filepath, 'rs') as fh:
        print("DADA File Header:")
        print("-" * 40)
        for key, value in fh.header.items():
            print(f"{key:20} {value}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path_to_dada_file>")
        sys.exit(1)
    
    dada_file = sys.argv[1]
    try:
        print_dada_header(dada_file)
    except Exception as e:
        print(f"Error reading DADA file: {e}")
