#include<iostream>
#include<vector>
#include<random>
#include<fstream>
#include<algorithm>
class Service{
public:
  static bool randomBit(double prob){
     static std::default_random_engine generator(std::random_device{}());
     static std::bernoulli_distribution distribution(prob); //0.5 really just means I want 50% true and 50% false, might be anything
      return distribution(generator);  
  }
  static int randomInt(int MAX_INT){
     std::mt19937 mt{}; //oldest trick in the book: take mod x to keep in [0. x)
     return mt()%MAX_INT;
  }
  static double f(std::vector<double> p, double x){
     return p[0]*x*x+p[1]*x+p[2];
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
        return c1.getFitness()<c2.getFitness();
     }
     void runIteration(){///the magic happens here
         ///STEP ONE: generate fitness scores for all chromosomes
         this->setChromosomes();
         for(auto c: genome)
           c.printBits();
        std::vector<Chromosome> newGenome;
         ///STEP TWO: select 1-crossover_ratio-mutate_ratio of the base set for progressing them, since they're clearly good
         int index;
         std::sort(genome.begin(), genome.end(),  compareChromosomes);
         for(index=0; index<(int)(pop_size*(1-crossover_ratio-mutate_ratio)); index++){
             
         }
         //these are THE SPECIAL ONEs, and they advance to the next change unmoved
         ///STEP THREE: crossover_ratio times, pick 2 of THE SPECIAL ONES and recombine them (THE KIDS)
         while(index<(int)(pop_size*(1-mutate_ratio))){
            //we choose who breeds
            int dad=Service::randomInt(index-1), mom=Service::randomInt(index-1);
            int breakingPoint=Service::randomInt(chromo_size);
            index++;
         }
         ///STEP FOUR: From among all these, flip a few positions on some and check

     }
};
int main(){
    std::ifstream fin("hyperpar.in");
    int pop_size=30, chromo_size=15, lo=0, hi=2, stages=1;
    double crossoverRatio=0.2, mutateRatio=0.1, a=1, b=-1, c=0;
    std::vector<double> params{a, b, c};
    Operation env(pop_size, chromo_size, lo, hi, stages, crossoverRatio, mutateRatio);
    env.setParameters(params);
    env.runIteration();
    return 0;
}