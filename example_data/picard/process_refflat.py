"""Download the current GENCODE data in refFlat format.

Connects to UCSC Genome Browser MySQL instance.

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A

"""


from sqlalchemy import *
import csv
import argparse

def main():
    db = create_engine('mysql://genome@genome-mysql.cse.ucsc.edu/hg19')
    db.echo = False  # Try changing this to True and see what happens
    metadata = MetaData(db)

    # The users table already exists, so no need to redefine it. Just
    # load it from the database using the "autoload" feature.
    gencode_table = Table('wgEncodeGencodeBasicV7', metadata, autoload=True)

    s = select([gencode_table.c.name, gencode_table.c.name2, gencode_table.c.chrom,
                gencode_table.c.strand, gencode_table.c.txStart, gencode_table.c.txEnd,
                gencode_table.c.cdsStart, gencode_table.c.cdsEnd, gencode_table.c.exonCount,
                gencode_table.c.exonStarts, gencode_table.c.exonEnds])
    
    rs = s.execute()
    
    of = file("/tmp/hg19_wgEncodeGencodeBasicV7_20120222.refFlat", "w")
    wr = csv.writer(of, delimiter="\t")

    try:
        while True:
            data = rs.fetchmany(1000)
            wr.writerows(data)
    except:
        of.close()

if __name__ == "__main__":
    main()

