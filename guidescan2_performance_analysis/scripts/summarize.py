import argparse
import glob
import re
import os.path

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Summarizes timing results into a CSV format."
    )

    parser.add_argument(
        "timing_dirs", help="Directories containing timing results for each tool",
        nargs='+'
    )

    return parser.parse_args()

def process_benchmark(tool_name, benchmark):
    filename = os.path.splitext(os.path.basename(benchmark))[0]

    m = re.search(
        r'(?P<guidetype>[^_]+)_(?P<sample_size>\d+)_k(?P<mismatches>\d+)_t(?P<threads>\d+)_(?P<pam>.*)',
        filename
    )

    guide_type = m.group('guidetype')
    sample_size = m.group('sample_size')
    num_mismatches = m.group('mismatches')
    threads = m.group('threads')
    pam = m.group('pam')

    with open(benchmark, 'r') as f:
        benchmark_results = dict(zip(next(f).split(), next(f).split()))

    seconds = benchmark_results['s']
    max_memory_mb = benchmark_results['max_rss']
    cpu_time = benchmark_results['cpu_time']

    return ','.join([tool_name, guide_type, sample_size, num_mismatches, threads, pam, seconds, cpu_time, max_memory_mb])

if __name__ == "__main__":
    args = parse_arguments()

    print('Tool,gRNA Set,Sample Size,Mismatches,Threads,PAM,Wall Clock Time (s),CPU Time (s),Max Memory (MB)')
    for timing_dir in args.timing_dirs:
        tool_name = timing_dir.split('/')[1]
        benchmarks = glob.glob(timing_dir + "/*.benchmark")
        csv_rows = [process_benchmark(tool_name, b) for b in benchmarks]
        print('\n'.join(csv_rows))
