#!/usr/bin/env python
import argparse
import sys

import cluster_commands
import try_command


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('job_ids',
                        nargs='+',
                        help='the cluster id(s) to cancel')
    args = parser.parse_args()
    return args


def run_cancel_command(command):
    retry_interval_seconds = list()  # no retries
    stdout, error = try_command.try_command(command, retry_interval_seconds)
    if error:
        sys.exit(error)

    return stdout


def main():
    args = parse_args()
    job_ids = args.job_ids
    command = cluster_commands.cancel_command(job_ids)
    run_cancel_command(command)


if __name__ == '__main__':
    main()
