#pragma once // shame on you! #pragma twice, shame on me!

#include <vector>

#include "defs/Track.h"

namespace analysis
{
	namespace core
	{
		class Flag : public Track
		{
			public:
				Flag() : Track() {this->reset();}

				virtual void reset()
				{

					_goodVertices = false;
					_globalSuperTightHalo2016Filter = false;
					_EcalDeadCellTriggerPrimitiveFilter = false;
					_BadPFMuonFilter = false;
					_BadPFMuonDzFilter = false;
					_hfNoisyHitsFilter = false;
					_eeBadScFilter = false;
					_ecalBadCalibFilter = false;
					Track::reset();
					_track.reset();

				}
				virtual ~Flag() {}

				bool _goodVertices;
				bool _globalSuperTightHalo2016Filter;
				bool _EcalDeadCellTriggerPrimitiveFilter;
				bool _BadPFMuonFilter;
				bool _BadPFMuonDzFilter;
				bool _hfNoisyHitsFilter;
				bool _eeBadScFilter;
				bool _ecalBadCalibFilter;
				Track _track;
		};
		typedef std::vector<analysis::core::Flag> Flags;
	}
}
