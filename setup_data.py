"""

This file exists to save genomes in order to analyze them in main.py.
Please run file and paste genome info from chosen database as per instructions.

"""

import sqlite3
import os

def write_genome(path, bases):
    # Connects to database
    connection = sqlite3.connect(path)
    cursor = connection.cursor()

    # Creates a table for genome
    cursor.execute("CREATE TABLE IF NOT EXISTS genome(base_sequence TEXT)")

    # Inserts bases
    cursor.execute("INSERT OR REPLACE INTO genome VALUES (?)", (bases,))

    # Passes data
    connection.commit()
    connection.close()

def create_genome():
    # Accesses dataset folder
    folder = "Genomes"
    os.makedirs(folder, exist_ok=True)  # <-- creates folder if it doesn't exist

    # Makes a .db for the organism to map
    name = input("What is the name of the organism:\n").replace(' ', '_')
    db_path = os.path.join(folder, f"{name}-genome.db")

    print(f'\nPlease paste bases')
    print("Press ENTER on empty line when complete.\n")

    # Loops through all lines on the pasted dataset
    lines = []
    while True:
        line = input()
        if line == '':
            break
        lines.append(line)

    # Joins all lines together
    bases = ''.join(lines)

    # Removes all invalid characters
    bases = ''.join(c.upper() for c in bases if c.lower() in {'a', 't', 'c', 'g', 'n'})

    write_genome(db_path, bases)

if __name__ == '__main__':
    create_genome()
    # TODO: ADD OPTIONS FOR MAPPING PROTEINS TO FURTHER DEBUG IN main.py