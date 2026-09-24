#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <map>
#include <iterator>
#include <cstdint>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include "AmQStrTdcData_RM.h"
#include "FileSinkHeader.h"
#include "FileSinkTrailer.h"
#include "FilterHeader.h"
#include "HeartbeatFrameHeader.h"
#include "SubTimeFrameHeader.h"
#include "TimeFrameHeader.h"

using namespace std;

enum {
	kFileSinkHeader,
  kTimeFrame,
  kSubTimeFrame,
  kFilter,
  kHeartbeatFrame,
  kFileSinkTrailer,
  kNumHeaderType
};
enum {
  kAmQData,
  kAmQHeartbeat,
  kAmQHeartbeat2nd,
  kAmQTrailer,
  kAmQThrottlingT1Start,
  kAmQThrottlingT1End,
  kAmQThrottlingT2Start,
  kAmQThrottlingT2End,
  kAmQStreamingRM
};
enum {
  kHRTDC = 5,
  kLRTDC = 6,
  kNumFETypes = 2
};

constexpr uint64_t kMagicWords[kNumHeaderType] = {FileSinkHeader::v1::MAGIC,
																									TimeFrame::v1::MAGIC,
                                                  SubTimeFrame::v1::MAGIC,
																									Filter::v1::MAGIC,
																									HeartbeatFrame::MAGIC,
																									FileSinkTrailer::v1::MAGIC};
constexpr int kDataSize[kNumHeaderType] = {38, 3, 6, 7, 2, 38};

constexpr uint32_t kHbFrameMask = 0x00ffffff;
constexpr uint32_t kHbFrameHalfRange = 0x00800000;

// Track observed intervals until EOF instead of diagnosing adjacent arrivals.
// Use signed, unwrapped positions: delayed frames may precede the first arrival.
// Assumption: reordering/discontinuities are less than half the 24-bit cycle.
struct FrameCoverage {
  bool initialized = false;
  int64_t highest = 0;
  std::map<int64_t, int64_t> ranges;
};
unordered_map<uint32_t, FrameCoverage> frameCoverageByFem;

void CheckFrameNumber(uint32_t femId, uint32_t frameNumber) {
  frameNumber &= kHbFrameMask;
  auto& coverage = frameCoverageByFem[femId];
  int64_t position = frameNumber;
  if (coverage.initialized) {
    const uint32_t delta =
        (frameNumber - (static_cast<uint32_t>(coverage.highest) & kHbFrameMask))
        & kHbFrameMask;
    position = coverage.highest + (delta < kHbFrameHalfRange
        ? static_cast<int64_t>(delta)
        : static_cast<int64_t>(delta) - (int64_t{1} << 24));
  } else {
    coverage.initialized = true;
    coverage.highest = position;
  }
  if (position > coverage.highest) coverage.highest = position;

  auto& ranges = coverage.ranges;
  auto next = ranges.upper_bound(position);
  auto current = ranges.end();
  if (next != ranges.begin()) {
    auto previous = std::prev(next);
    if (position <= previous->second) return; // duplicate META frame
    if (position == previous->second + 1) {
      previous->second = position;
      current = previous;
    }
  }
  if (current == ranges.end())
    current = ranges.emplace(position, position).first;
  if (next != ranges.end() && next->first == current->second + 1) {
    current->second = next->second;
    ranges.erase(next);
  }
}

void ReportFrameGaps() {
  // Only internal gaps in the observed input are knowable. Do not infer missing
  // frames before the first or after the last observed frame for an FE.
  for (const auto& entry : frameCoverageByFem) {
    const auto& ranges = entry.second.ranges;
    if (ranges.empty()) continue;
    auto previous = ranges.begin();
    for (auto current = std::next(previous); current != ranges.end(); ++current) {
      const auto delta = current->first - previous->second;
      cerr << Form(
          "Frame gap (EOF coverage): FEM=%08x previous=%u current=%u",
          entry.first,
          static_cast<uint32_t>(previous->second) & kHbFrameMask,
          static_cast<uint32_t>(current->first) & kHbFrameMask)
          << " delta=" << delta << " missing=" << (delta - 1) << endl;
      previous = current;
    }
  }
}

void PrintThrottlingEvent(
    int eventType,
    uint32_t femId,
    uint32_t timeFrameId,
    const char* section,
    int wordIndex) {
  const char* eventName = nullptr;
  switch (eventType) {
    case kAmQThrottlingT1Start: eventName = "T1 start"; break;
    case kAmQThrottlingT1End:   eventName = "T1 end";   break;
    case kAmQThrottlingT2Start: eventName = "T2 start"; break;
    case kAmQThrottlingT2End:   eventName = "T2 end";   break;
    default: return;
  }

  cout << Form(
      "Throttling %s: section=%s FEM=%08x IP=%u.%u.%u.%u TF=%u word=%d",
      eventName,
      section,
      femId,
      (femId >> 24) & 0xff,
      (femId >> 16) & 0xff,
      (femId >> 8) & 0xff,
      femId & 0xff,
      timeFrameId,
      wordIndex)
      << endl;
}

bool read_u64_le(istream& is, uint64_t& out) {
  unsigned char b[8];
  if (!is.read(reinterpret_cast<char*>(b), 8)) return false;
  out = (uint64_t)b[0]
     | ((uint64_t)b[1] << 8)
     | ((uint64_t)b[2] << 16)
     | ((uint64_t)b[3] << 24)
     | ((uint64_t)b[4] << 32)
     | ((uint64_t)b[5] << 40)
     | ((uint64_t)b[6] << 48)
     | ((uint64_t)b[7] << 56);
  return true;
}

int CheckHeader(uint64_t buf) {
	for (int i = 0; i < kNumHeaderType; ++i) {
		if (buf == kMagicWords[i]) {
			return i;
		}
	}
	return -1;
}

void ReadExact(istream* ifs, void* dst, size_t n) {
  if (!ifs->read(reinterpret_cast<char*>(dst), (streamsize)n))
    throw runtime_error("EOF/IO error");
}

int CheckAMANEQHeader(AmQStrTdc::Data::v1::Bits bits) {
  namespace AmQ = AmQStrTdc::Data::v1;
  auto head = bits.head;
  if (head == AmQ::HeadTypes::Data) {
    return kAmQData;
  } else if (head == AmQ::HeadTypes::Heartbeat) {
    return kAmQHeartbeat;
  } else if (head == AmQ::HeadTypes::Heartbeat2nd) {
    return kAmQHeartbeat2nd;
  } else if (head == AmQ::HeadTypes::Trailer) {
    return kAmQTrailer;
  } else if (head == AmQ::HeadTypes::ThrottlingT1Start) {
    return kAmQThrottlingT1Start;
  } else if (head == AmQ::HeadTypes::ThrottlingT1End) {
    return kAmQThrottlingT1End;
  } else if (head == AmQ::HeadTypes::ThrottlingT2Start) {
    return kAmQThrottlingT2Start;
  } else if (head == AmQ::HeadTypes::ThrottlingT2End) {
    return kAmQThrottlingT2End;
  } else if (head == AmQ::HeadTypes::StreamingRM) {
    return kAmQStreamingRM;
  }
  return -1;
}

int ReadMeta(int n_fe, ifstream* ifs) {
  uint64_t buf;
  SubTimeFrame::v1::Header  h_stf{};
  HeartbeatFrame::Header    h_hbf{};
	AmQStrTdc::Data::v1::Bits h_amq{};
  for (int i = 0; i < n_fe; ++i) {
    // SUBTIME
    read_u64_le(*ifs, buf);
	  auto flag = CheckHeader(buf);
    if (flag != kSubTimeFrame) {
      cerr << Form("Unexpected flag (ReadMeta): %d != %d", flag, kSubTimeFrame) << endl;
      return -1;
    }
    ReadExact(ifs, reinterpret_cast<char*>(&h_stf) + sizeof(uint64_t), sizeof(SubTimeFrame::v1::Header) - sizeof(uint64_t));
    // HRTBEAT
    read_u64_le(*ifs, buf);
	  flag = CheckHeader(buf);
    if (flag != kHeartbeatFrame) {
      cerr << Form("Unexpected flag (ReadMeta): %d != %d", flag, kHeartbeatFrame) << endl;
      return -1;
    }
    ReadExact(ifs, reinterpret_cast<char*>(&h_hbf) + sizeof(uint64_t), sizeof(HeartbeatFrame::Header) - sizeof(uint64_t));
    auto length = h_hbf.length / 8 - kDataSize[kHeartbeatFrame];
    for (int k = 0; k < length; ++k) {
      read_u64_le(*ifs, buf);
      memcpy(&h_amq, &buf, 8);
      auto amq_flag = CheckAMANEQHeader(h_amq);
      PrintThrottlingEvent(
          amq_flag, h_stf.femId, h_stf.timeFrameId, "META", k);
      // A throttled STF can contain additional heartbeat words in its META
      // payload.  Only the first word of the final heartbeat delimiter pair
      // represents the HBF being summarized here; do not mix embedded frames
      // into the continuity check.
      const bool isFinalHeartbeat =
          amq_flag == kAmQHeartbeat && k + 2 == length;
      if (isFinalHeartbeat) {
        CheckFrameNumber(h_stf.femId, h_amq.hbframe);
      }
    }
  }
  return 0;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    cerr << "./bin [input.dat] [output.root]" << endl;
    return -1;
  }
  auto in_name  = TString(argv[1]);
  auto out_name = TString(argv[2]);

  ifstream ifs(in_name);
  if (!ifs) {
    cerr << "failed to open " << in_name << endl;
    return -1;
  }

	auto out_file = TFile(out_name, "recreate");
	auto out_tree = TTree("tree", "tree");
// filter
  uint16_t flt_type;
  uint32_t n_triggers;
  uint32_t worker_id;
  uint32_t elapse_time;
  uint64_t flt_sec;
  uint64_t flt_usec;
  uint16_t tth_type;
  uint32_t trg_time;
  vector<uint32_t> trg_type;
// hit
  int n_hits;
  vector<uint32_t> fem_type;
  vector<uint32_t> ip_addr;
  vector<uint8_t>  ip2;
  vector<uint8_t>  ip3;
//  vector<uint32_t> frame_number;
//  vector<uint32_t> laccp_offset;
//  vector<uint32_t> hb_flags;
  uint32_t frame_number;
  uint32_t laccp_offset;
  uint32_t hb_flags;
  vector<uint64_t> tdc;
  vector<uint64_t> tot;
  vector<uint64_t> ch;
  vector<uint64_t> data_type;
  out_tree.Branch("flt_type",     &flt_type,    "flt_type/s");
  out_tree.Branch("n_triggers",   &n_triggers,  "n_triggers/i");
  out_tree.Branch("worker_id",    &worker_id,   "worker_id/i");
  out_tree.Branch("elapse_time",  &elapse_time, "elapse_time/i");
  out_tree.Branch("flt_sec",      &flt_sec,     "flt_sec/l");
  out_tree.Branch("flt_usec",     &flt_usec,    "flt_usec/l");
  out_tree.Branch("tth_type",     &tth_type,    "tth_type/s");
  out_tree.Branch("trg_time",     &trg_time,    "trg_time/i");
  out_tree.Branch("trg_type",     &trg_type);
  out_tree.Branch("fem_type",     &fem_type);
//  out_tree.Branch("frame_number", &frame_number);
//  out_tree.Branch("laccp_offset", &laccp_offset);
//  out_tree.Branch("hb_flags",     &hb_flags);
  out_tree.Branch("frame_number", &frame_number, "frame_number/i");
  out_tree.Branch("laccp_offset", &laccp_offset, "laccp_offset/i");
  out_tree.Branch("hb_flags",     &hb_flags,     "hb_flags/i");
  out_tree.Branch("n_hits",       &n_hits,       "n_hits/I");
  out_tree.Branch("ip_addr",      &ip_addr);
  out_tree.Branch("ip2",          &ip2);
  out_tree.Branch("ip3",          &ip3);
  out_tree.Branch("tdc",          &tdc);
  out_tree.Branch("tot",          &tot);
  out_tree.Branch("ch",           &ch);
  out_tree.Branch("data_type",    &data_type);
  
  uint64_t buf;
  uint32_t tf_id_from_tf = 0;
  FileSinkHeader::v1::Header   h_fsh{};
  TimeFrame::v1::Header        h_tf{};
  SubTimeFrame::v1::Header     h_stf{};
	Filter::v1::Header           h_flt{};
  Filter::v1::TrgTimeHeader    h_tth{};
  Filter::v1::TrgTime          h_tt{};
  HeartbeatFrame::Header       h_hbf{};
	FileSinkTrailer::v1::Trailer h_fst{};
	AmQStrTdc::Data::v1::Bits    h_amq{};
  int tmp_ip_addr;
  int tmp_ip2;
  int tmp_ip3;
  vector<uint32_t> tmp_trg_time;
  int n_tf = 0;
  while (read_u64_le(ifs, buf)) {
		auto flag = CheckHeader(buf);
//cerr << "data: " << hex << buf << endl;
//cerr << "flag: " << dec << flag << endl;
		if (flag == -1) {
			cerr << "Invalid magic word: " << hex << buf << endl;
			break;
		}
		if (flag == kFileSinkHeader) {
    	ReadExact(&ifs, reinterpret_cast<char*>(&h_fsh) + sizeof(uint64_t), sizeof(FileSinkHeader::v1::Header) - sizeof(uint64_t));
		  cout << "FairMQ device type: " << h_fsh.fairMQDeviceType << endl;
		  cout << "Run number: "         << h_fsh.runNumber        << endl;
		  tm tm = *localtime(&h_fsh.startUnixtime);
      cout << "Start: " << put_time(&tm, "%Y-%m-%d %H:%M:%S %Z") << endl;
    } else if (flag == kTimeFrame) {
    	ReadExact(&ifs, reinterpret_cast<char*>(&h_tf) + sizeof(uint64_t), sizeof(TimeFrame::v1::Header) - sizeof(uint64_t));
      auto n_fe = h_tf.numSource;
      // FLTLOGI
      read_u64_le(ifs, buf);
		  flag = CheckHeader(buf);
      if (flag != kFilter) {
        cerr << Form("Unexpected flag: %d != %d, position: %d", flag, kFilter, ifs.tellg()) << endl;
        return -1;
      }
      ReadExact(&ifs, reinterpret_cast<char*>(&h_flt) + sizeof(uint64_t), sizeof(Filter::Header) - sizeof(uint64_t));
      flt_type    = h_flt.type;
//           tf_id       = h_flt.timeFrameId;
      n_triggers  = h_flt.numTrigs;
      worker_id   = h_flt.workerId;
//      auto n_messages  = h_flt.numMessages;
      elapse_time = h_flt.elapseTime;
      flt_sec     = h_flt.timeSec;
      flt_usec    = h_flt.timeUSec;
      tmp_trg_time.clear();
      trg_type.clear();
      tmp_trg_time.resize(n_triggers);
      trg_type.resize(n_triggers);
//cerr << Form("%d %d %d %d %d %d %d %d %d %d", h_flt.length, h_flt.hLength, h_flt.type, h_flt.timeFrameId, h_flt.numTrigs, h_flt.workerId, h_flt.numMessages, h_flt.elapseTime, h_flt.timeSec, h_flt.timeUSec) << endl;
      // TRGTIME
      ReadExact(&ifs, reinterpret_cast<char*>(&h_tth), sizeof(Filter::TrgTimeHeader));
      tth_type = h_tth.type;
//cerr << Form("%d %d %d", h_tth.length, h_tth.hLength, h_tth.type) << endl;
      for (int i = 0; i < n_triggers; ++i) {
        read_u64_le(ifs, buf);
        memcpy(&h_tt, &buf, 8);
        tmp_trg_time[i] = h_tt.time * 4;
        trg_type[i] = h_tt.type & 0x00ffffffu;
			}
      ReadMeta(n_fe, &ifs);
      for (int i = 0; i < n_triggers; ++i) {
        fem_type.clear();
        ip_addr.clear();
        ip2.clear();
        ip3.clear();
        tdc.clear();
        tot.clear();
        ch.clear();
        data_type.clear();
        // TIMEFRM
        read_u64_le(ifs, buf);
		    flag = CheckHeader(buf);
        if (flag != kTimeFrame) {
          cerr << Form("Unexpected flag: %d != %d", flag, kTimeFrame) << endl;
          return -1;
        }
    	  ReadExact(&ifs, reinterpret_cast<char*>(&h_tf) + sizeof(uint64_t), sizeof(TimeFrame::v1::Header) - sizeof(uint64_t));
        n_fe = h_tf.numSource;
        for (int j = 0; j < n_fe; ++j) {
          // SUBTIME
          read_u64_le(ifs, buf);
		      flag = CheckHeader(buf);
          if (flag != kSubTimeFrame) {
            cerr << Form("Unexpected flag: %d != %d", flag, kSubTimeFrame) << endl;
            return -1;
          }
    	    ReadExact(&ifs, reinterpret_cast<char*>(&h_stf) + sizeof(uint64_t), sizeof(SubTimeFrame::v1::Header) - sizeof(uint64_t));
          fem_type.push_back(h_stf.femType);
          tmp_ip_addr = h_stf.femId;
          tmp_ip2 = (tmp_ip_addr >> 8) & 0xFF;
          tmp_ip3 = (tmp_ip_addr >> 0) & 0xFF;
          // HRTBEAT
          read_u64_le(ifs, buf);
		      flag = CheckHeader(buf);
          if (flag != kHeartbeatFrame) {
            cerr << Form("Unexpected flag: %d != %d", flag, kHeartbeatFrame) << endl;
            return -1;
          }
          ReadExact(&ifs, reinterpret_cast<char*>(&h_hbf) + sizeof(uint64_t), sizeof(HeartbeatFrame::Header) - sizeof(uint64_t));
          auto length = h_hbf.length / 8 - kDataSize[kHeartbeatFrame];
          for (int k = 0; k < length; ++k) {
            read_u64_le(ifs, buf);
            memcpy(&h_amq, &buf, 8);
            auto amq_flag = CheckAMANEQHeader(h_amq);
            if (amq_flag == kAmQData) {
              if (fem_type.back() == kLRTDC) {
                tdc.push_back(h_amq.tdc);
                tot.push_back(h_amq.tot);
                ch.push_back(h_amq.ch);
                data_type.push_back(h_amq.dtype);
              } else if (fem_type.back() == kHRTDC) {
                tdc.push_back(h_amq.hrtdc);
                tot.push_back(h_amq.hrtot);
                ch.push_back(h_amq.hrch);
                data_type.push_back(h_amq.hrdtype);
              }
              ip_addr.push_back(tmp_ip_addr);
              ip2.push_back(tmp_ip2);
              ip3.push_back(tmp_ip3);
            } else if (amq_flag == kAmQHeartbeat) {
              hb_flags     = h_amq.hbflag;
              laccp_offset = h_amq.toffset;
              frame_number = h_amq.hbframe;
            } else if (amq_flag == kAmQHeartbeat2nd) {
            } else if (amq_flag == kAmQThrottlingT1Start ||
                       amq_flag == kAmQThrottlingT1End ||
                       amq_flag == kAmQThrottlingT2Start ||
                       amq_flag == kAmQThrottlingT2End) {
              PrintThrottlingEvent(
                  amq_flag, tmp_ip_addr, h_stf.timeFrameId, "SLICE", k);
            } else if (amq_flag == -1) {
              cerr << "Invalid header: " << hex << h_amq.head << endl;
              return -1;
            }
          }
        }
        n_hits = tdc.size();
        trg_time = tmp_trg_time[i];
        out_tree.Fill();
      }
		} else if (flag == kFileSinkTrailer) {
    	ReadExact(&ifs, reinterpret_cast<char*>(&h_fst) + sizeof(uint64_t),
	  					  sizeof(FileSinkTrailer::v1::Trailer) - sizeof(uint64_t));
		  tm tm = *localtime(&h_fst.stopUnixtime);
      cout << "Stop: " << put_time(&tm, "%Y-%m-%d %H:%M:%S %Z") << endl;
		}
	}
  ReportFrameGaps();
  out_file.Write();
  return 0;
}
