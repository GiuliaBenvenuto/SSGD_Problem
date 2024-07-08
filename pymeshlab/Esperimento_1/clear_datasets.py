import csv
import re
import sys

def is_very_small_number(s):
    # Check if the string represents a very small number (e.g., 2.37358e-313)
    pattern = r'^-?\d+\.?\d*e-\d{3,}$'
    return bool(re.match(pattern, s))

def process_csv(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        for row in reader:
            new_row = []
            for cell in row:
                if is_very_small_number(cell.strip()):
                    new_row.append('0')
                else:
                    new_row.append(cell)
            writer.writerow(new_row)

    print(f"Processed {input_file} and saved results to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    process_csv(input_file, output_file)