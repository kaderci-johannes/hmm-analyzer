#pragma once

#include <vector>

#include "defs/Track.h"

namespace analysis
{
	namespace core
	{
		class Photon : public Track
		{
			public:
				Photon() : Track() {this->reset();}

				virtual void reset()
				{
					_hoe = 0;
					_isScEtaEB;
					_isScEtaEE;
					_mvaID;
					_pixelSeed;
					Track::reset();
					_track.reset();

				}
				virtual ~Photon() {}

				Track _track;
				float _hoe;
				bool _isScEtaEB;
				bool _isScEtaEE;
				float _mvaID;
				bool _pixelSeed;


		};

		typedef std::vector<analysis::core::Photon> Photons;
	}
}
