#include<iostream>
#include<vector>
#include<random>
#include<fstream>
#include<algorithm>
#include<map>
#include<string>
class Service{
public:
  static std::ofstream fout;
  static bool randomBit(double prob){
     static std::default_random_engine generator(std::random_device{}());
     static std::bernoulli_distribution distribution(prob); //0.5 really just means I want 50% true and 50% false, might be anything
      return distribution(generator);  
  }
  static double modul(double x){
    return (x<0) ? -x : x;
  }
  static double randomUniform(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distro(0.0, 1.0);
    return distro(gen);
  }
  static int randomInt(int MAX_INT){
     std::mt19937 mt{}; //oldest trick in the book: take mod x to keep in [0. x)
     return mt()%MAX_INT;
  }
  /// @brief  
  // This is the function we're optimizing. For now, no negative values in range, it fucks with the averaging 
  static double f(std::vector<double> p, double x){
     return p[0]*x*x+p[1]*x+p[2];
  } //TODO: find some way to map all real values into [0, 1] 
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
std::ofstream Service::fout("log.txt");
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
        ///eg. magic numbers, 
        for(bool g: genes)
            Service::fout<<g;
        Service::fout<<"\n";
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
    bool showMechanism;
    std::vector<double> param; 
    public:
     Operation(int ps, int cs, int lo_, int hi_, int stages_, double crossr, double mutr):
     pop_size(ps), chromo_size(cs), lo(lo_), hi(hi_), stages(stages_), crossover_ratio(crossr), mutate_ratio(mutr){
         showMechanism=true;
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
     void runIteration(int i){
    /// Printing initial genome if i == 0
    std::string outputFile("log.txt");
    if(i == 0){
        Service::fout << "Populatia initiala:\n";
        for(auto &x : genome)
            x.printBits();
    }

    /// STEP ONE: generate fitness scores for all chromosomes
    std::vector<Chromosome> newGenome;
    double sum = 0;
    std::vector<double> relativeFitness, absoluteFitness;
    for(auto &chromosome : genome){
        absoluteFitness.push_back(sum);
        sum += chromosome.getFitness(); //This cannot take negative vlaues, they don;t easily map to [0, 1]
    }
    for(double fitScore : absoluteFitness)
            relativeFitness.push_back(fitScore / sum);
    relativeFitness.push_back(1.0);  //flag for ease of computation
    /// Print partial sums and probability of selection if i == 0
    if(i == 0){
        Service::fout << "Probabilitate de selectie per cromozom:\n";
        Service::fout<<relativeFitness[0]<<"\n";
        for(size_t j = 1; j < relativeFitness.size(); j++){
            Service::fout<<relativeFitness[j]-relativeFitness[j-1]<<"\n";
        }
        Service::fout<< "Probablitate cumulativa: \n";
        for(size_t j=0; j<relativeFitness.size()-1; j++){
           Service::fout<<relativeFitness[j]<<"\n";
        }
    }

    /// STEP TWO: Select chromosomes based on fitness
    int index;
    std::sort(genome.begin(), genome.end(), compareChromosomes);
    newGenome.push_back(genome[0]);

    for(index = 1; index < (int)(pop_size * (1 - crossover_ratio - mutate_ratio)); index++){
        double seed = Service::randomUniform();
        int targetChromo = Service::binarySearch(relativeFitness, seed);
        if(i==0){
          Service::fout<<"Cu seedul "<<seed;
          Service::fout<<" alegem sa progresam cromozomul "<<targetChromo<<"\n";
        }
        newGenome.push_back(genome[targetChromo]);
    }///we picked a bunch of THE SPECIAL ONES, these go forward as is

    /// Print genome after selection if i == 0
    if(i == 0){
        Service::fout << "Dupa doua etape\n";
        for(auto &x : newGenome)
            x.printBits();
    }

    /// STEP THREE: In meiosis, homologues cross over each other and informtion changes spots
    while(index < (int)(pop_size * (1 - mutate_ratio))){
        int dad = Service::randomInt(index - 1);
        int mom = Service::randomInt(index - 1);
        int breakingPoint = Service::randomInt(chromo_size);

        std::vector<bool> new1, new2;
        for(int j = 0; j < chromo_size; j++){
            if(j < breakingPoint){
                new1.push_back(newGenome[dad].getGene(j)); //the way this basically works is 
                new2.push_back(newGenome[mom].getGene(j));
            } else {
                new1.push_back(newGenome[mom].getGene(j));
                new2.push_back(newGenome[dad].getGene(j));
            }
        }

        if(i == 0){
            Service::fout << "Am facut incrucisarea lui ";
            newGenome[dad].printBits();
            Service::fout << " cu ";
            newGenome[mom].printBits(); 
            Service::fout<<" si am taiat la " << breakingPoint << "\n";

        }

        newGenome.push_back(Chromosome(new1));
        newGenome.push_back(Chromosome(new2));
        Service::fout<<"Si se adauga: \n";
        newGenome[newGenome.size()-2].printBits();
        newGenome[newGenome.size()-1].printBits();
        index += 2;
    }
    ///We also output the genome after step 3, mutations
    if(i == 0){
        Service::fout << "Dupa pasul 3\n";
        for(auto &x : newGenome)
            x.printBits();
    }
    /// STEP FOUR: Mutation operation
    while(index < pop_size){
        int base = Service::randomInt(index - 1);
        Chromosome c1(genome[base]);
        c1.radiate();
        newGenome.push_back(c1);
        index++;
    }

    /// Print new genome if i == 0
    if(i == 0){
        Service::fout << "Si asta e populatia finala\n";
        for(auto &x : newGenome)
            x.printBits();
    }
    /// STEP FIVE: Overwrite the whole generation
    std::swap(this->genome, newGenome);
   }
};
int main(){
    std::ifstream fin("hyperpar.txt");
    int pop_size, chromo_size, lo, hi, stages;
    double crossoverRatio, mutateRatio, a, b, c;
    //it's ugly as fuck, but we have to read a lot of stuff with specific names
    fin>>pop_size>>chromo_size>>lo>>hi>>stages>>crossoverRatio>>mutateRatio>>a>>b>>c;
    std::vector<double> params{a, b, c};
    Operation env(pop_size, chromo_size, lo, hi, stages, crossoverRatio, mutateRatio);
    env.setParameters(params);
    env.setChromosomes();///we generate an initial population
    for(int i=0; i<stages; i++){
       env.runIteration(i);
       Service::fout<<"Generatia "<<i+1<<"maxim: "<<env.getMaxFitness()<<"\n";
       Service::fout<<"Generatia "<<i+1<<"medie: "<<env.getAverageFitness()<<"\n";
    }
    return 0;
}