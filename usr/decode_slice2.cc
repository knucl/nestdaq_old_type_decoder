#include <fstream>
#include <iomanip>
#include <iostream>

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

int pastFrameNumber = -1;

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
      cerr << "Unexpected flag: " << flag << endl;
      return -1;
    }
    ReadExact(ifs, reinterpret_cast<char*>(&h_stf) + sizeof(uint64_t), sizeof(SubTimeFrame::v1::Header) - sizeof(uint64_t));
    // HRTBEAT
    read_u64_le(*ifs, buf);
	  flag = CheckHeader(buf);
    if (flag != kHeartbeatFrame) {
      cerr << "Unexpected flag: " << flag << endl;
      return -1;
    }
    ReadExact(ifs, reinterpret_cast<char*>(&h_hbf) + sizeof(uint64_t), sizeof(HeartbeatFrame::Header) - sizeof(uint64_t));
    auto length = h_hbf.length / 8 - kDataSize[kHeartbeatFrame];
    for (int k = 0; k < length; ++k) {
      read_u64_le(*ifs, buf);
      memcpy(&h_amq, &buf, 8);
      auto amq_flag = CheckAMANEQHeader(h_amq);
      if (amq_flag == kAmQHeartbeat) {
        auto frame_number = h_amq.hbframe;
        if (pastFrameNumber != -1 &&
            frame_number != pastFrameNumber &&
            frame_number != pastFrameNumber + 1) {
          cerr << Form("Invalid frame number: %d (expect: %d or +1)", frame_number, pastFrameNumber) << endl;
        }
        pastFrameNumber = frame_number;
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
//  vector<uint32_t> trg_type;
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
//  out_tree.Branch("trg_type",     &trg_type);
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
        cerr << "Unexpected flag: " << flag << endl;
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
//      trg_type.clear();
      tmp_trg_time.resize(n_triggers);
//      trg_type.resize(n_triggers);
//cerr << Form("%d %d %d %d %d %d %d %d %d %d", h_flt.length, h_flt.hLength, h_flt.type, h_flt.timeFrameId, h_flt.numTrigs, h_flt.workerId, h_flt.numMessages, h_flt.elapseTime, h_flt.timeSec, h_flt.timeUSec) << endl;
      // TRGTIME
      ReadExact(&ifs, reinterpret_cast<char*>(&h_tth), sizeof(Filter::TrgTimeHeader));
      tth_type = h_tth.type;
//cerr << Form("%d %d %d", h_tth.length, h_tth.hLength, h_tth.type) << endl;
      for (int i = 0; i < n_triggers; ++i) {
        read_u64_le(ifs, buf);
        memcpy(&h_tt, &buf, 8);
        tmp_trg_time[i] = h_tt.time * 4;
//        trg_type[i] = h_tt.type;
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
          cerr << "Unexpected flag: " << flag << endl;
          return -1;
        }
    	  ReadExact(&ifs, reinterpret_cast<char*>(&h_tf) + sizeof(uint64_t), sizeof(TimeFrame::v1::Header) - sizeof(uint64_t));
        n_fe = h_tf.numSource;
        for (int j = 0; j < n_fe; ++j) {
          // SUBTIME
          read_u64_le(ifs, buf);
		      flag = CheckHeader(buf);
          if (flag != kSubTimeFrame) {
            cerr << "Unexpected flag: " << flag << endl;
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
            cerr << "Unexpected flag: " << flag << endl;
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
            } else if (amq_flag == kAmQThrottlingT1Start) {
              cout << "Throttling T1 start" << endl;
              cout << "IP: " << int{tmp_ip3} << endl;
            } else if (amq_flag == kAmQThrottlingT1End) {
              cout << "Throttling T1 end"   << endl;
              cout << "IP: " << int{tmp_ip3} << endl;
            } else if (amq_flag == kAmQThrottlingT2Start) {
              cout << "Throttling T2 start" << endl;
              cout << "IP: " << int{tmp_ip3} << endl;
            } else if (amq_flag == kAmQThrottlingT2End) {
              cout << "Throttling T2 end"   << endl;
              cout << "IP: " << int{tmp_ip3} << endl;
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
  out_file.Write();
  return 0;
}
