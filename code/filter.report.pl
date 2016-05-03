$drugname=$ARGV[0];
$num=$ARGV[1];
open OUTPUT, ">gsea_report_${drugname}.txt" or die $!;
print OUTPUT "Gene\tES\tNES\tNom p-val\tFDR q-val\tFWER p-val\n";

open FILE, "${drugname}.Gsea.${num}/gsea_report_for_vehicle_${num}.xls" or die $!;
<FILE>;
while(<FILE>){
	chomp;
	$line1=$_;
	@line=split("\t",$line1);
	print OUTPUT "$line[1]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\n";
}
close FILE;

open FILE, "${drugname}.Gsea.${num}/gsea_report_for_drug_${num}.xls" or die $!;
<FILE>;
while(<FILE>){
	chomp;
	$line1=$_;
	@line=split("\t",$line1);
	print OUTPUT "$line[1]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\n";
}
close FILE;

close OUTPUT;
