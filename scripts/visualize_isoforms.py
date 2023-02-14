import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Visualize structure and abundance of isoforms'))
    parser.add_argument(
        '--gene-id',
        required=True,
        help='The gene_id to visualize')
    parser.add_argument('--gtf',
                        required=True,
                        help='The path to the ESPRESSO output updated.gtf')
    parser.add_argument('--out-dir',
                        required=True,
                        help='The path to use as output directory')

    return parser.parse_args()


def visualize_isoforms(args):
    # TODO
    raise Exception('visualize_isoforms: not implemented')


def main():
    args = parse_args()
    visualize_isoforms(args)


if __name__ == '__main__':
    main()
