import sys


def trim_quotes(string):
    return string.strip('"')


def read_genes_from_gtf(path):
    by_id = dict()
    by_name = dict()
    with open(path, 'rt') as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('#'):
                continue

            columns = line.split('\t')
            attribute_string = columns[8]
            attrs = attribute_string.split(';')
            gene_id = ''
            gene_name = ''
            for attr in attrs:
                attr = attr.strip()
                if not attr:
                    continue

                parts = attr.split()
                if parts[0] == 'gene_id':
                    gene_id = trim_quotes(parts[1])
                if parts[0] == 'gene_name':
                    gene_name = trim_quotes(parts[1])

            if gene_id != '':
                by_id[gene_id] = gene_name
            if gene_name != '':
                by_name[gene_name] = gene_id

    return {'by_id': by_id, 'by_name': by_name}


def find_gene_id(gene_name, gene_id, updated_gtf, ref_gtf):
    id_and_name_maps_updated = read_genes_from_gtf(updated_gtf)
    id_and_name_maps_ref = read_genes_from_gtf(ref_gtf)
    if gene_id != '':
        found = id_and_name_maps_updated['by_id'].get(gene_id)
        if found is not None:
            return gene_id

        found = id_and_name_maps_ref['by_id'].get(gene_id)
        if found is not None:
            return gene_id

    if gene_name != '':
        found = id_and_name_maps_updated['by_name'].get(gene_name, '')
        if found != '':
            return found

        found = id_and_name_maps_ref['by_name'].get(gene_name, '')
        if found != '':
            return found

    return None


def main():
    gene_name = sys.argv[1]
    gene_id = sys.argv[2]
    updated_gtf = sys.argv[3]
    ref_gtf = sys.argv[4]
    found = find_gene_id(gene_name, gene_id, updated_gtf, ref_gtf)
    if found:
        print(found)
    else:
        print('')


if __name__ == '__main__':
    main()
