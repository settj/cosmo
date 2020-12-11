#include <iostream>
#include <fstream>

void getKmers( size_t& nKmers, size_t k, std::set< std::string >& kmers, std::string& p) {
   std::ifstream in;
   in.open(p);
   if(in.fail())
    {
        std::cout<<"Could not read the fasta file called "<<p<<std::endl;
        exit(0);
    }
   std::string sline;
   std::vector< std::string > vline;
   if (!getline(in, sline) || !(sline[0] == '>')) {
          std::cout << "Make sure that "<<p<< " is a fatsa file" << std::endl;
          exit(0);
        }
   while ( getline( in, sline ) ) {
      vline.push_back( sline );
   }

   size_t pos = 0;
   std::vector< std::string > reads;
   std::string read;
   do {
      if (vline[pos][0] == '>') {
	 //finish current read and start a new one
	 if (!read.empty()) {
	    reads.push_back(read);
	    read.clear();
	 }
      } else {
	 read += vline[pos];
      }
      ++pos;
   } while (pos != vline.size());

   if (!read.empty()) //handle the last read
      reads.push_back( read );

   for (size_t i = 0; i < reads.size(); ++i) {
      std::string sline = reads[i];
      size_t read_length = sline.size();
      if (read_length < k){
        std::cout<<"Warning! There is a read in file "<<p<<" that is shorter than k="<<k<<"!"<<std::endl;
        exit(0);
      }

      for (size_t i = 0; i < read_length; ++i) {
	       if (sline[i] == 'N'){
	        sline[i] = 'A';
        }
      }
      size_t nMers = read_length - k + 1;
      for (size_t start = 0; start < nMers; ++start) {
	       std::string kmer = sline.substr( start, k );
	           kmers.insert( kmer );
      }
   }

   in.close();
   nKmers = kmers.size();
}
