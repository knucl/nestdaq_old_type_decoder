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

constexpr int kDataHit      = 0;
// 1 - 5 are used in streaming RM
constexpr int kDataSpillOn  = 6;
constexpr int kDataSpillEnd = 7;

constexpr int kNumHBF = 4;

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
	uint64_t event_id = 0xffffffffffffffff;
  uint64_t spill_id = 0xffffffffffffffff;
  uint32_t tf_id;
	uint32_t hbf_id;
  uint32_t fem_type;
  uint32_t ip_addr;
  vector<uint8_t> ip(4);
//  uint32_t frame_number = 0;
  uint32_t laccp_offset = 0;
  uint64_t tdc;
  uint64_t tot;
  uint64_t ch;
  uint64_t data_type;
  out_tree.Branch("event_id",         &event_id,         "event_id/l");
  out_tree.Branch("spill_id",         &spill_id,         "spill_id/l");
  out_tree.Branch("tf_id",            &tf_id,            "tf_id/i");
  out_tree.Branch("hbf_id",           &hbf_id,           "hbf_id/i");
  out_tree.Branch("fem_type",         &fem_type,         "fem_type/i");
//  out_tree.Branch("frame_number", &frame_number, "frame_number/i");
  out_tree.Branch("laccp_offset",     &laccp_offset,     "laccp_offset/i");
  out_tree.Branch("ip_addr",          &ip_addr,          "ip_addr/i");
  out_tree.Branch("ip",               &ip);
  out_tree.Branch("tdc",              &tdc,              "tdc/l");
  out_tree.Branch("tot",              &tot,              "tot/l");
  out_tree.Branch("ch",               &ch,               "ch/l");
  out_tree.Branch("data_type",        &data_type,        "data_type/l");

  auto out_tree1 = TTree("tree1", "tree1");
  uint32_t hbf_id_flt;
  int n_triggers;
  vector<uint32_t> trg_time;
  out_tree1.Branch("tf_id",      &tf_id,      "tf_id/i");
  out_tree1.Branch("hbf_id",     &hbf_id_flt, "hbf_id/i");
  out_tree1.Branch("n_triggers", &n_triggers, "n_triggers/I");
  out_tree1.Branch("trg_time",   &trg_time);
  
  auto out_tree2 = TTree("tree2", "tree2"); // heart beat delimiter
  uint64_t hbd_flag;
  uint64_t trans_size;
  uint64_t gen_size;
  out_tree2.Branch("tf_id",      &tf_id,      "tf_id/i");
  out_tree2.Branch("hbf_id",     &hbf_id,     "hbf_id/i");
  out_tree2.Branch("ip",         &ip);
  out_tree2.Branch("hbd_flag",   &hbd_flag,   "hbd_flag/l");
  out_tree2.Branch("trans_size", &trans_size, "trans_size/l");
  out_tree2.Branch("gen_size",   &gen_size,   "gen_size/l");

//  auto out_tree3 = TTree("tree3", "tree3"); // can merge to tree (tree0)
//  int throttling_id;
//  uint32_t start_tf_id;
//  uint32_t start_hbf_id;
//  uint32_t end_tf_id;
//  uint32_t end_hbf_id;
//  out_tree3.Branch("ip",            &ip);
//  out_tree3.Branch("throttling_id", &throttling_id, "throttling_id/I");
//  out_tree3.Branch("start_tf_id",   &start_tf_id,   "start_tf_id/i");
//  out_tree3.Branch("start_hbf_id",  &start_hbf_id,  "start_hbf_id/i");
//  out_tree3.Branch("end_tf_id",     &end_tf_id,     "end_tf_id/i");
//  out_tree3.Branch("end_hbf_id",    &end_hbf_id,    "end_hbf_id/i");

  uint64_t buf;
  uint32_t length;
  int n_hbf = 0;
  bool have_tf_header = false;
  uint32_t tf_id_from_tf = 0;
  map<uint32_t, int> is_throttling; // IP, (T2, T1)
  FileSinkHeader::v1::Header   h_fsh{};
  TimeFrame::v1::Header        h_tf{};
  SubTimeFrame::v1::Header     h_stf{};
  Filter::v1::Header           h_flt{};
  Filter::v1::TrgTimeHeader    h_tth{};
  Filter::v1::TrgTime          h_tt{};
	HeartbeatFrame::Header       h_hbf{};
	FileSinkTrailer::v1::Trailer h_fst{};
	AmQStrTdc::Data::v1::Bits    h_amq{};
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
      tf_id_from_tf  = h_tf.timeFrameId;
      tf_id          = h_tf.timeFrameId;
      have_tf_header = true;
      n_hbf          = 0;
      hbf_id         = 0;
//h_tf.Print();
		} else if (flag == kSubTimeFrame) {
    	ReadExact(&ifs, reinterpret_cast<char*>(&h_stf) + sizeof(uint64_t), sizeof(SubTimeFrame::v1::Header) - sizeof(uint64_t));
      if (have_tf_header && tf_id_from_tf != h_stf.timeFrameId) {
        cerr << "TimeFrame ID mismatch: TF=" << tf_id_from_tf << " STF=" << h_stf.timeFrameId << endl;
      }
      tf_id    = h_stf.timeFrameId;
      fem_type = h_stf.femType;
      ip_addr  = h_stf.femId;
      ip[0] = (ip_addr >> 24) & 0xFF;
      ip[1] = (ip_addr >> 16) & 0xFF;
      ip[2] = (ip_addr >>  8) & 0xFF;
      ip[3] = (ip_addr >>  0) & 0xFF;
      n_hbf  = 0;
      hbf_id = 0;
//cerr << Form("FEM type: %u, IP addr.: ", h_stf.femType)
//     << ((ip_addr >> 24) & 0xFF) << '.'
//     << ((ip_addr >> 16) & 0xFF) << '.'
//     << ((ip_addr >>  8) & 0xFF) << '.'
//     << ((ip_addr >>  0) & 0xFF) << endl;
    } else if (flag == kFilter) {
//      if (!is_first) {
//        out_tree.Fill();
//      }
//      is_first = false;
      ReadExact(&ifs, reinterpret_cast<char*>(&h_flt) + sizeof(uint64_t), sizeof(Filter::Header) - sizeof(uint64_t));
//      flt_type    = h_flt.type;
////           tf_id       = h_flt.timeFrameId;
//      n_triggers  = h_flt.numTrigs;
//      worker_id   = h_flt.workerId;
////      auto n_messages  = h_flt.numMessages;
//      elapse_time = h_flt.elapseTime;
//      flt_sec     = h_flt.timeSec;
//      flt_usec    = h_flt.timeUSec;
//      trg_time.clear();
////      trg_type.clear();
//      trg_time.resize(n_triggers);
//      trg_type.resize(n_triggers);
//cerr << Form("%d %d %d %d %d %d %d %d %d %d", h_flt.length, h_flt.hLength, h_flt.type, h_flt.timeFrameId, h_flt.numTrigs, h_flt.workerId, h_flt.numMessages, h_flt.elapseTime, h_flt.timeSec, h_flt.timeUSec) << endl;
      for (int i = 0; i < kNumHBF; ++i) {
        hbf_id_flt = i;
        ReadExact(&ifs, reinterpret_cast<char*>(&h_tth), sizeof(Filter::TrgTimeHeader));
        n_triggers = h_tth.length / 8 - 2;
        trg_time.clear();
        trg_time.resize(n_triggers);
//cerr << Form("%d %d %d", h_tth.length, h_tth.hLength, h_tth.type) << endl;
        for (int j = 0; j < n_triggers; ++j) {
          read_u64_le(ifs, buf);
          memcpy(&h_tt, &buf, 8);
          trg_time[j] = h_tt.time;
//cerr << Form("%lu %lu", h_tt.time, h_tt.type) << endl;
//        trg_time[i] = h_tt.time * 4;
//        trg_type[i] = h_tt.type;
        }
        out_tree1.Fill();
      }
//      fem_type.clear();
//      ip_addr.clear();
//      ip2.clear();
//      ip3.clear();
//      frame_number.clear();
//      laccp_offset.clear();
//      hb_flags.clear();
//      tdc.clear();
//      tot.clear();
//      ch.clear();
//      data_type.clear();
    } else if (flag == kHeartbeatFrame) {
      ReadExact(&ifs, reinterpret_cast<char*>(&h_hbf) + sizeof(uint64_t), sizeof(HeartbeatFrame::Header) - sizeof(uint64_t));
      
      hbf_id = n_hbf;
      ++n_hbf;
//cerr << h_hbf.length << endl;
      length = h_hbf.length / 8 - kDataSize[kHeartbeatFrame];
//cerr << "length: " << length << endl;
      for (int i = 0; i < length; ++i) {
        read_u64_le(ifs, buf);
        memcpy(&h_amq, &buf, 8);
//cerr << hex << buf << endl;
        auto amq_flag = CheckAMANEQHeader(h_amq);
//cerr << amq_flag << endl;
        if (amq_flag == kAmQData) {
          if (fem_type == kLRTDC) {
            tdc       = h_amq.tdc;
            tot       = h_amq.tot;
            ch        = h_amq.ch;
            data_type = kDataHit;
          } else if (fem_type == kHRTDC) {
            tdc       = h_amq.hrtdc;
            tot       = h_amq.hrtot;
            ch        = h_amq.hrch;
            data_type = kDataHit;
          } else {
            cerr << "Inlvaid FE type" << endl;
            continue;
          }
          out_tree.Fill();
        } else if (amq_flag == kAmQHeartbeat) {
          hbd_flag     = h_amq.hbflag;
          laccp_offset = h_amq.toffset;
//          frame_number = h_amq.hbframe;
        } else if (amq_flag == kAmQHeartbeat2nd) {
          trans_size = h_amq.transSize;
          gen_size   = h_amq.geneSize;
          if (hbf_id == 3) {
            out_tree2.Fill();
          }
        } else if (amq_flag == kAmQThrottlingT1Start) {
          cout << "Throttling T1 start" << endl;
          cout << "IP: " << int{ip[3]} << endl;
        } else if (amq_flag == kAmQThrottlingT1End) {
          cout << "Throttling T1 end"   << endl;
          cout << "IP: " << int{ip[3]} << endl;
        } else if (amq_flag == kAmQThrottlingT2Start) {
          cout << "Throttling T2 start" << endl;
          cout << "IP: " << int{ip[3]} << endl;
        } else if (amq_flag == kAmQThrottlingT2End) {
          cout << "Throttling T2 end"   << endl;
          cout << "IP: " << int{ip[3]} << endl;
        } else if (amq_flag == kAmQStreamingRM) {
          tdc       = h_amq.rmtdc;
          tot       = 0xffffffffffffffff;
          ch        = 0xffffffffffffffff;
          event_id  = h_amq.eventID;
          spill_id  = h_amq.spillID;
          data_type = h_amq.rmch;
          out_tree.Fill();
        } else {
          cerr << "Invalid Amq Header: " << hex << h_amq.head << endl;
          return -1;
        }
      }
		} else if (flag == kFileSinkTrailer) {
    	ReadExact(&ifs, reinterpret_cast<char*>(&h_fst) + sizeof(uint64_t),
	  					  sizeof(FileSinkTrailer::v1::Trailer) - sizeof(uint64_t));
		  tm tm = *localtime(&h_fst.stopUnixtime);
      cout << "Stop: " << put_time(&tm, "%Y-%m-%d %H:%M:%S %Z") << endl;
		} else {
      cerr << "Invalid header" << endl;
      return -1;
    }
	}

  out_file.Write();
  return 0;
}
