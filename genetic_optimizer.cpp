#include<iostream>
#include<vector>
#include<random>
#include<fstream>
#include<algorithm>
#include<map>
class Service{
public:
  static bool randomBit(double prob){
     static std::default_random_engine generator(std::random_device{}());
     static std::bernoulli_distribution distribution(prob); //0.5 really just means I want 50% true and 50% false, might be anything
      return distribution(generator);  
  }
  static bool randomUniform(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distro(0.0, 1.0);
    return distro(gen);
  }
  static int randomInt(int MAX_INT){
     std::mt19937 mt{}; //oldest trick in the book: take mod x to keep in [0. x)
     return mt()%MAX_INT;
  }
  static double f(std::vector<double> p, double x){
     return p[0]*x*x+p[1]*x+p[2];
  }
  static int binarySearch(std::vector<double> &v, double target){
        int step=(1<<16), n=v.size(), currentPos=0;
        while(step>0){
            if(currentPos+step<n && v[currentPos+step]<target)
               currentPos+=step;
            step>>=1;
        }
        return currentPos;
    }
};
class Chromosome{
    /// @brief "environmental" stuff: parameters for objective function, and interval over which we optimize
    /// which is commmon to all chromosomes
    static double lo, hi;
    static std::vector<double> params; ///maybe there are more than 3 parameters
    std::vector<bool> genes;
public:
    Chromosome(int n){
        for(int  i=0; i<n; i++){
            genes.push_back(Service::randomBit(0.5));
        }
    }
    Chromosome(const std::vector<bool> &g): genes(g) {}
    //we use 'ionizing radiation' on the chromosome, to make a new one(we can't destroy the elder)
    void radiate(){ 
        int position=Service::randomInt(genes.size());
        genes[position]=(!genes[position]);
    }
    bool getGene(int locus){
        return genes[locus];
    }
    static void setEnvironment(double lo_, double hi_, std::vector<double>&& param_){
        lo=lo_, hi=hi_;
        params=std::move(param_);
    }
    //Turns the genes of the chromosome into the actual point represented 
    double decode(){
        double value=lo, step=(hi-lo)/2;
        for(auto gene: genes){
            if(gene==1)
               value+=step;
            step/=2;  
        }
        return value;
    }
    ///this just applies the polynomial (or any other function, for that matter) on the actual value
    double getFitness(){
       return Service::f(this->params, this->decode()); 
    }
    void printBits(){
        for(bool g: genes)
            std::cout<<g;
        std::cout<<"\n";
    }
};
//sucks to do this, but I must initalize the env
double Chromosome::lo=0.0, Chromosome::hi=0.0;
std::vector<double> Chromosome::params={};
/// @brief 
// A class with a single method, which runs the steps in the actual algorithm
class Operation{  
    std::vector<Chromosome> genome;
    int pop_size, chromo_size, lo, hi, stages;
    double crossover_ratio, mutate_ratio;
    std::vector<double> param; 
    public:
     Operation(int ps, int cs, int lo_, int hi_, int stages_, double crossr, double mutr):
     pop_size(ps), chromo_size(cs), lo(lo_), hi(hi_), stages(stages_), crossover_ratio(crossr), mutate_ratio(mutr) {

     }
     void setParameters(std::vector<double> &ps){
        for(auto parameter: ps)
          param.push_back(parameter);
     }
     void setChromosomes(){
         for(int i=0; i<pop_size; i++){
            genome.emplace_back(Chromosome(chromo_size));
         }
         Chromosome::setEnvironment(lo, hi, std::move(param));
     }
     static double compareChromosomes(Chromosome c1, Chromosome c2){
        return c1.getFitness()>c2.getFitness();
     }
     double getMaxFitness(){
        double currentMaxFitness=-1e7;
        for(auto chrom: genome)
          currentMaxFitness=std::max(currentMaxFitness, chrom.getFitness());
        return currentMaxFitness;
     }
     double getAverageFitness(){
        double sumFitness=0;
        for(auto chrom: genome)
           sumFitness+=chrom.getFitness();
        return sumFitness/pop_size;
     }
     void runIteration(int i){///the magic happens here
         ///STEP ONE: generate fitness scores for all chromosomes
        std::vector<Chromosome> newGenome;
        for(auto x: genome)
            x.printBits();
        ///STEP 1.2: build a method of selecting chromosomes proportional to fitness
        double sum=0;
        std::vector<double> relativeFitness, absoluteFitness;
        for(auto chromosome: genome){
            absoluteFitness.push_back(sum);
            sum+=chromosome.getFitness();
        }
        
        for(double fitScore: absoluteFitness){
            relativeFitness.push_back(fitScore/sum); ///TODO see if we can avoid dividing floats, those take way tooo long
            ///TRICK: the ith entry in the array is always the i-th chromosome
        } 
         ///STEP TWO: select 1-crossover_ratio-mutate_ratio of the base set for progressing them, since they're clearly good
         int index;
         std::sort(genome.begin(), genome.end(),  compareChromosomes);
         newGenome.push_back(genome[0]);
         //we will do so via "roulette selection", pick a random number in the array and see what we get
         for(index=1; index<(int)(pop_size*(1-crossover_ratio-mutate_ratio)); index++){
             double seed=Service::randomUniform();
             int targetChromo=Service::binarySearch(relativeFitness, seed);
             newGenome.push_back(genome[targetChromo]);
         }
         //these are THE SPECIAL ONEs, and they advance to the next change unmoved
         ///STEP THREE: crossover_ratio times, pick 2 of THE SPECIAL ONES and recombine them (THE KIDS)
         while(index<(int)(pop_size*(1-mutate_ratio))){
            //we choose who breeds
            int dad=Service::randomInt(index-1), mom=Service::randomInt(index-1);
            int breakingPoint=Service::randomInt(chromo_size);
            std::vector<bool> new1, new2;
            for(int i=0; i<chromo_size; i++){
                if(i<breakingPoint){
                   new1.push_back(newGenome[dad].getGene(i));
                   new2.push_back(newGenome[mom].getGene(i));
                }
                else{
                   new1.push_back(newGenome[mom].getGene(i));
                   new2.push_back(newGenome[dad].getGene(i));
                }
            }
            newGenome.push_back(Chromosome(new1));
            newGenome.push_back(Chromosome(new2));
            index+=2;
         }
         ///STEP FOUR: From among all these, flip a few positions on some and check
        while(index<pop_size){
            int base=Service::randomInt(index-1);
            Chromosome c1(genome[base]);
            c1.radiate();
            newGenome.push_back(c1);
            index++;///FINALLY, IT LEARNS
        }
        ///STEP FIVE: overwrite the whole generation
        std::swap(this->genome, newGenome);
     }
};
int main(){
    std::ifstream fin("hyperpar.in");
    int pop_size=30, chromo_size=4, lo=-5, hi=15, stages=2;
    double crossoverRatio=0.2, mutateRatio=0.1, a=1, b=-1, c=0;
    std::vector<double> params{a, b, c};
    Operation env(pop_size, chromo_size, lo, hi, stages, crossoverRatio, mutateRatio);
    env.setParameters(params);
    env.setChromosomes();///we generate an initial population
    for(int i=0; i<stages; i++){
       env.runIteration(i);
       std::cout<<"Generatia "<<i+1<<"maxim: "<<env.getMaxFitness()<<"\n";
       std::cout<<"Generatia "<<i+1<<"medie: "<<env.getAverageFitness()<<"\n";
    }
    return 0;
}