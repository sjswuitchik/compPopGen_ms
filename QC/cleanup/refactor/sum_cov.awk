#!/usr/bin/awk -f

# usage when not called in write_coverage_beds.sh: gzip -dc spp_name.merge.bg.gz | awk -v avg=x -v spp=spp_name -f sum_cov.awk

BEGIN { 
    FS = OFS = "\t"
    print "chrom", "start", "end", "cov" > spp"_coverage_sites_low.bed"
    print "chrom", "start", "end", "cov" > spp"_coverage_sites_high.bed"
    print "chrom", "start", "end", "cov" > spp"_coverage_sites_clean.bed"
}
{
    mean = avg
    if ($4 < 0.5*mean)
        print $1, $2, $3, $4 > spp"_coverage_sites_low.bed"
    else if ($4 > 2.0*mean)
        print $1, $2, $3, $4 > spp"_coverage_sites_high.bed"
    else
        print $1, $2, $3, $4 > spp"_coverage_sites_clean.bed"
}
