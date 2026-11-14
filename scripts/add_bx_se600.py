import pysam
import sys

def main(input_bam, output_bam):
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
            for read in infile:
                if read.has_tag("BX"):
                    outfile.write(read)
                    continue

                # readname: @...#1_2_3/1  or  @...#1_2_3
                try:
                    bc = read.query_name.split("#")[1].split("/")[0]
                    if bc != "0_0_0":  # 0_0_0 = no valid barcode; skip to avoid fake linked-fragment in extractHAIRS
                        read.set_tag("BX", bc, value_type="Z")
                except (IndexError, AttributeError):
                    pass

                outfile.write(read)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} input.bam output.bam")
    main(sys.argv[1], sys.argv[2])