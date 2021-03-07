kitties=nancycats()
dataset=kitties.loci
#a=allele_freq(testlocus.genotype[3])
#b=allele_freq(testlocus.genotype[4])
#popA=testlocus[testlocus[!,:population].=="1",:]
#popgenotype=[]
totalalleles = unique(genlocus[!,:locus])
SS=DataFrame(popid=String[],ssdpop=Float32[])

for g in totalalleles
    genelocus=data[data[!,:locus].==g,:]
    ssd=Euclidean(dropmissing(genelocus))
    ssdtotal=ssd/(2*length(eachrow(genelocus)))
    push!(SS,["ssd(t)",sstotal])

    rp=unique(genelocus[!,:population]) # holds the names of all the populations
    for k in rp #interate through all populations
        popgenes=dropmissing(genelocus[genelocus[!, :population].==k,:]) # extract population k from the whole dataset
        ss= Euclidean(popgenes)#calculate the ssd for k population
        push!(SS,[k ss])
    end
end



function Euclidean(popgenotype::DataFrame)
    deltasum=0
    popallelecount=length(eachrow(popgenotype)) # number of individuals in population k
    #@time begin
    for i in 1:popallelecount-1 #iterate through all individuals in population k -1
        j=i+1 #counter for genotype of the next individual
        a=allele_freq(popgenotype.genotype[i])# get the allele frequency of the first individual 
        while j<=popallelecount #check to see if we have interated through all of the other individuals in the population
            b=allele_freq(popgenotype.genotype[j]) # get the allele frequencies of the second individual
            difference=values(mergewith(-,a,b))# subtract the allele frequencies between the two individuals
            deltasum=deltasum+sum(x->x^2,difference) # square the differences
            j+=1 # add one to the second individual counter
        end
    end
return deltasum
end
