#pragma once

#include <vector>

#include "defs/Track.h"

namespace analysis
{
	namespace core
	{
		class TrigObj : public Track
		{
			public:
				TrigObj() : Track() {this->reset();}

				virtual void reset()
				{
					Track::reset();
					_track.reset();

					_id = 0;
					_filterBits = 0;
				}
				virtual ~TrigObj() {}

				Track _track;

				ushort _id;
				int _filterBits;

		};

		typedef std::vector<analysis::core::TrigObj> TrigObjs;
	}
}
