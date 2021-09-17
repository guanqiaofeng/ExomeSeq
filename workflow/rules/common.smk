
intervals = pd.read_table(
    config["bed_file"]
).set_index(
    "interval", drop=False
)
def get_intervals(wildcards):
    inter = wildcards.interval
    bed = "/cluster/home/selghamr/workflows/ExomeSeq/resources/hg38_bed/" + inter + ".bed"
    return bed


def get_MuTect2_output(wildcards):
    outfile = "{sample}_{interval}.mut2.vcf"
    res = []
        res.append(
            "results/MuTect2/{}/{}_{}".format(
                sample, interval, outfile
            )
        )
    return res
