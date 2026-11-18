#ifndef Analysis_Core_Track_h
#define Analysis_Core_Track_h

#include "Object.h"

namespace analysis
{
	namespace core
	{
		class Track : public Object
		{
			public:
				Track() : Object() {this->reset();}

				virtual void reset()
				{
					_charge = 0;
					_pt = 0;
					_pterr = 0;
					_eta = 0;
					_phi = 0;
					_mass = 0;
					_ptcorr = 0;
					_etacorr = 0;
					_phicorr = 0;
					_masscorr = 0;
				}
				virtual ~Track() {}

				int _charge;
				float _pt;
				float _pterr;
				float _eta;
				float _phi;
				float _mass;
				float _ptcorr;
				float _etacorr;
				float _phicorr;
				float _masscorr;

		};

		typedef std::vector<analysis::core::Track> Tracks;
	}
}

#endif
