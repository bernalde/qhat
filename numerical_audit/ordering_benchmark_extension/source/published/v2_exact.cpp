// Independent implementation of Tranter Eq. (11), coefficient-l1 variant.
// Input: number of terms, then x z integer_coefficient dense_lex_label.
// Input terms arrive in descending magnitude, ties by dense lex label.
// argv[1]: append-only checkpoint, one accepted insertion per line.
#include <boost/multiprecision/cpp_int.hpp>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <unistd.h>
using boost::multiprecision::cpp_int;
using I=cpp_int;
struct Term {uint32_t x,z; I c;int label,y;};
inline int pc(uint32_t x){return __builtin_popcount(x);}
inline bool anti(uint32_t x,uint32_t z,uint32_t a,uint32_t b){return pc((x&b)^(z&a))&1;}
inline int phase(const Term&a,const Term&b){return (a.y+b.y-pc((a.x^b.x)&(a.z^b.z))+2*pc(a.z&b.x))&3;}
inline I ab(const I&a){return a<0?-a:a;}
struct Poly {std::unordered_map<uint64_t,I> v;I norm=0;
 void add(uint64_t k,const I&d){auto it=v.find(k);if(it==v.end()){if(d!=0){v.emplace(k,d);norm+=ab(d);}return;}norm-=ab(it->second);it->second+=d;if(it->second==0)v.erase(it);else norm+=ab(it->second);}
};
void outer(Poly&p,const Term&a,uint32_t x,uint32_t z,int innerph,const I&innercoef,int factor){
 if(!anti(a.x,a.z,x,z))return;
 int ph=(innerph+a.y+pc(x&z)-pc((a.x^x)&(a.z^z))+2*pc(a.z&x))&3;
 if(ph&1)throw std::runtime_error("Nonreal double commutator");
 I val=a.c*innercoef;val*=factor;if(ph==2)val=-val;
 p.add((uint64_t)(a.x^x)|((uint64_t)(a.z^z)<<32),val);
}
int main(int argc,char**argv){
 try {
 if(argc!=2)throw std::runtime_error("checkpoint path required");int n;std::cin>>n;std::vector<Term> ts(n);for(auto&t:ts){std::cin>>t.x>>t.z>>t.c>>t.label;t.y=pc(t.x&t.z);}if(!std::cin)throw std::runtime_error("Invalid input");
 std::vector<int> order;std::ifstream prior(argv[1]);int label,pos;std::string norm;
 while(prior>>label>>pos>>norm){int next=order.size();if(next>=n||label!=ts[next].label||pos<0||pos>next)throw std::runtime_error("Checkpoint does not match input");order.insert(order.begin()+pos,next);}
 if(!prior.eof()&&prior.fail())throw std::runtime_error("Partial or corrupt checkpoint");
 FILE*out=fopen(argv[1],"a");if(!out)throw std::runtime_error("Cannot open checkpoint");
 for(int ni=order.size();ni<n;++ni){const Term&N=ts[ni];Poly p;p.v.reserve(100000);
  // Appending N: delta W=[2P+N,[N,P]], W=-24 V2.
  for(int ci:order){const Term&C=ts[ci];if(!anti(N.x,N.z,C.x,C.z))continue;
   auto x=N.x^C.x,z=N.z^C.z;int ph=phase(N,C);I co=4*N.c*C.c;
   for(int ai:order)outer(p,ts[ai],x,z,ph,co,2);outer(p,N,x,z,ph,co,1);
  }
  I best=p.norm;int bestpos=order.size();
  for(int at=(int)order.size()-1;at>=0;--at){const Term&A=ts[order[at]];
   if(anti(N.x,N.z,A.x,A.z)){auto x=N.x^A.x,z=N.z^A.z;int ph=phase(N,A);I co=4*N.c*A.c;
    // Moving N left across A: -6[P,[N,A]]-3[A+N,[N,A]].
    for(int k=0;k<at;++k)outer(p,ts[order[k]],x,z,ph,co,-6);
    outer(p,A,x,z,ph,co,-3);outer(p,N,x,z,ph,co,-3);
   }
   if(p.norm<=best){best=p.norm;bestpos=at;}
  }
  order.insert(order.begin()+bestpos,ni);std::string s=best.convert_to<std::string>();
  fprintf(out,"%d %d %s\n",N.label,bestpos,s.c_str());fflush(out);fsync(fileno(out));
  if((ni+1)%64==0||ni+1==n)std::cerr<<"accepted "<<ni+1<<"/"<<n<<"\n";
 }
 fclose(out);std::cout<<"[";bool first=true;for(int i:order){if(!first)std::cout<<",";first=false;std::cout<<ts[i].label;}std::cout<<"]\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
}
