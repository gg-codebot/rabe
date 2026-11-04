"""Write a copy of an electric_rabe restart file that looks like a killed run.

Keeps the header and all fieldlines before the last one, then half of the
lines of the last fieldline, and finally half of the following line.
"""

import sys


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <complete.restart.dat> <truncated.restart.dat>")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        lines = f.readlines()

    header = [line for line in lines if line.startswith("#")]
    data = [line for line in lines if not line.startswith("#")]
    ids = [int(line.split()[0]) for line in data]
    last_id = max(ids)
    if last_id < 2:
        print("FAIL: need at least two fieldlines in the restart file")
        sys.exit(1)

    kept = [line for line, i in zip(data, ids) if i < last_id]
    last = [line for line, i in zip(data, ids) if i == last_id]
    n_half = len(last) // 2
    cut_line = last[n_half][: len(last[n_half]) // 2]

    with open(sys.argv[2], "w") as f:
        f.writelines(header + kept + last[:n_half])
        f.write(cut_line)

    print(f"kept fieldlines 1..{last_id - 1}, cut fieldline {last_id} "
          f"after {n_half} of {len(last)} lines")


if __name__ == "__main__":
    main()
