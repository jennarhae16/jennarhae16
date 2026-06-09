import subprocess
import pandas as pd
import plotly.express as px

plink = 'plink2 --vcf assignment4.vcf --allow-extra-chr 0 --max-alleles 2 --make-bed --out assignment4'

subprocess.run(plink, shell = True)

ad2 = 'admixture -a qn2 --cv assignment4.bed 2 > log2.out'
ad3 = 'admixture -a qn2 --cv assignment4.bed 3 > log3.out'
ad4 = 'admixture -a qn2 --cv assignment4.bed 4 > log4.out'
ad5 = 'admixture -a qn2 --cv assignment4.bed 5 > log5.out'
ad6 = 'admixture -a qn2 --cv assignment4.bed 6 > log6.out'
ad7 = 'admixture -a qn2 --cv assignment4.bed 7 > log7.out'
ad8 = 'admixture -a qn2 --cv assignment4.bed 8 > log8.out'
ad9 = 'admixture -a qn2 --cv assignment4.bed 9 > log9.out'
ad10 = 'admixture -a qn2 --cv assignment4.bed 10 > log10.out'
ad11 = 'admixture -a qn2 --cv assignment4.bed 11 > log11.out'
ad_list = [ad2,ad3,ad4,ad5,ad6,ad7,ad8,ad9,ad10,ad11]
for i in ad_list:
    subprocess.run(i, shell = True)

cv_err1 = 'grep "CV" *out | awk \'{print $3,$4}\' | sed -e \'s/(//;s/)//;s/://;s/K=//\' > assignment4.cv.error'
cv_err2 = 'grep "CV" *out | awk \'{print $3,$4}\' | cut -c 4,7-20 > assignment4.cv.error'
cv_err3 = 'awk \'/CV/ {print $3,$4}\' *out | cut -c 4,7-20 > assignment4.cv.error'
cv_err4 = 'cat assignment4.cv.error'

cv_list = [cv_err1,cv_err2,cv_err3,cv_err4]

for j in cv_list:
    subprocess.run(j, shell= True)


assign4_data = input("Enter best K value file:\n")

df = pd.read_table(assign4_data, sep = " ", header = None)

headers = ["Pop1","Pop2","Pop3","Pop4","Pop5","Pop6","Pop7","Pop8","Pop9","Pop10","Pop11"]
df.columns = headers

df["Strain"]= df.index

fig = px.bar(df, x = "Strain", y = ["Pop1","Pop2","Pop3","Pop4","Pop5","Pop6","Pop7","Pop8","Pop9","Pop10","Pop11"], title = "Strains by Population", labels={"value": "Value", "variable": "Population"})
fig.show()


