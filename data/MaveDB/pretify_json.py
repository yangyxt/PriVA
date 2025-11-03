#!/usr/bin/env python
import json
import os

def prettify_json(input_path, output_path, indent=2):
    """
    Read a JSON file and write it back with proper indentation and formatting.
    
    Args:
        input_path: Path to the input JSON file
        output_path: Path where the formatted JSON will be saved
        indent: Number of spaces for indentation (default: 2)
    """
    print(f"Reading JSON from {input_path}...")
    
    try:
        with open(input_path, 'r') as f:
            data = json.load(f)
        
        print(f"Formatting JSON with indent={indent}...")
        
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_path)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=indent)
        
        print(f"Formatted JSON saved to {output_path}")
        print(f"Original file size: {os.path.getsize(input_path) / (1024*1024):.2f} MB")
        print(f"Formatted file size: {os.path.getsize(output_path) / (1024*1024):.2f} MB")
        
        # Print a small sample to show the formatting
        print("\nSample of the formatted JSON:")
        sample_size = 20  # Number of lines to show
        
        with open(output_path, 'r') as f:
            sample = []
            for i, line in enumerate(f):
                if i < sample_size:
                    sample.append(line.rstrip())
                else:
                    break
            
            print("\n".join(sample))
            print("...")
            
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    input_path = "data/MaveDB/meta.json"
    output_path = "data/MaveDB/meta_formatted.json"
    prettify_json(input_path, output_path) 