"""

This file exists to save genomes in order to analyze them in main.py.
Please run file and paste genome info from chosen database as per instructions.

"""

import sqlite3
import os

if __name__ == '__main__':
    # Accesses dataset folder
    folder = "Genomes"
    os.makedirs(folder, exist_ok=True)  # <-- creates folder if it doesn't exist

    # Makes a .db for the organism to map
    name = input("What is the name of the organism:\n").replace(' ', '_')
    db_path = os.path.join(folder, f"{name}-genome.db")

    # Gets total chromosome number
    num = 0

    ready = False
    while not ready:
        num = input('What is the chromosome number? ')
        try:
            # Checks for valid input
            num = int(num)
            ready = True
        except ValueError:
            print("Please enter valid integer.\n")

    # Gets chromosome values
    chromosomes = []

    for i in range(num):
        print(f'\nPlease paste bases in chromosome {i + 1}.')
        print("Press ENTER on empty line when complete.\n")

        # Loops through all lines on the pasted dataset
        lines = []
        while True:
            line = input()
            if line == '':
                break
            lines.append(line)

        # Joins all lines together
        chromosome = ''.join(lines)

        # Removes all invalid characters
        chromosome = ''.join(c.upper() for c in chromosome if c.lower() in {'a','t','c','g', 'n'})
        chromosomes.append(chromosome)

    # Connects to database
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()

    # Creates a table for genome
    cursor.execute("CREATE TABLE IF NOT EXISTS genome(chromosome_number INTEGER PRIMARY KEY, base_sequence TEXT)")

    # Inserts bases
    for i, seq in enumerate(chromosomes):
        cursor.execute("INSERT OR REPLACE INTO genome VALUES (?, ?)", (i + 1, seq))

    # Passes data
    connection.commit()
    connection.close()