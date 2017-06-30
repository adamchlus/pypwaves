"""Module :mod:`pypwaves.base` basic functions for opening and displaying pulsewave files"""


import struct, numpy as np,os, inspect
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from rtree import index




class PulseWaves(object):
    """Pulsewaves class object"""
    
    def __init__(self,pls_file):
        pulsebinary =  open(pls_file, 'rb')
              
        #read header
        self.filename = pls_file
        self.file_sig = "".join(struct.unpack("c"*16, pulsebinary.read(16))).strip("\x00")
        self.global_params =	struct.unpack("=L", pulsebinary.read(4))[0]
        self.file_id = struct.unpack("=L", pulsebinary.read(4))[0]
        self.proj_GUID1 = struct.unpack("=L", pulsebinary.read(4))[0]
        self.proj_GUID2 = struct.unpack("H", pulsebinary.read(2))[0]
        self.proj_GUID3 = struct.unpack("H", pulsebinary.read(2))[0]
        self.proj_GUID3 = struct.unpack("B"*8, pulsebinary.read(8))[0]
        self.sys_id ="".join(struct.unpack("c"*64, pulsebinary.read(64))).strip("\x00")
        self.software =	"".join(struct.unpack("c"*64, pulsebinary.read(64))).strip("\x00")
        self.file_day = struct.unpack("H", pulsebinary.read(2))[0]
        self.file_year=struct.unpack("H", pulsebinary.read(2))[0]
        self.version_maj = struct.unpack("B", pulsebinary.read(1))[0]
        self.version_min = struct.unpack("B", pulsebinary.read(1))[0]
        self.header_size = struct.unpack("H", pulsebinary.read(2))[0]
        self.offset_to_pulses = struct.unpack("q", pulsebinary.read(8))[0]
        self.num_pulses = struct.unpack("q", pulsebinary.read(8))[0]
        self.pulse_format = struct.unpack("=L", pulsebinary.read(4))[0]
        self.pulse_attr = struct.unpack("=L", pulsebinary.read(4))[0]
        self.pulse_size = struct.unpack("=L", pulsebinary.read(4))[0]
        self.pulse_compression =	struct.unpack("=L", pulsebinary.read(4))[0]
        self.reserved =	struct.unpack("q", pulsebinary.read(8))[0]
        self.num_vlr =	struct.unpack("I", pulsebinary.read(4))[0]
        self.num_avlr =	struct.unpack("!l", pulsebinary.read(4))[0]
        self.t_scale  =	struct.unpack("d", pulsebinary.read(8))[0]
        self.t_offset =	struct.unpack("d", pulsebinary.read(8))[0]
        self.t_min = struct.unpack("q", pulsebinary.read(8))[0]
        self.t_max = struct.unpack("q", pulsebinary.read(8))[0]
        self.x_scale = struct.unpack("d", pulsebinary.read(8))[0]
        self.y_scale = struct.unpack("d", pulsebinary.read(8))[0]
        self.z_scale = struct.unpack("d", pulsebinary.read(8))[0]
        self.x_offset =	struct.unpack("d", pulsebinary.read(8))[0]
        self.y_offset =	struct.unpack("d", pulsebinary.read(8))[0]
        self.z_offset =	struct.unpack("d", pulsebinary.read(8))[0]
        self.x_min = struct.unpack("d", pulsebinary.read(8))[0]
        self.x_max = struct.unpack("d", pulsebinary.read(8))[0]
        self.y_min = struct.unpack("d", pulsebinary.read(8))[0]
        self.y_max = struct.unpack("d", pulsebinary.read(8))[0]
        self.z_min = struct.unpack("d", pulsebinary.read(8))[0]
        self.z_max = struct.unpack("d", pulsebinary.read(8))[0]
        
        self.vlrs = {}
        self.avlrs = {}

        #read variable length records (VLR)
        for num_vlr in range(self.num_vlr):
            
            vlr = VLR(pulsebinary)

            #if vlr is a scanner
            if  vlr.record_id >=100001 and vlr.record_id < 100255:    
                vlr.record=Scanner(pulsebinary)
        
            #if vlr is a pulse descriptor 
            elif vlr.record_id >=200001 and vlr.record_id < 200255:              
                #read pulse desciptor
                vlr.record = PulseDecriptor(pulsebinary)            
                vlr.sampling_records = {}                
                #read sampling record
                for x in range(vlr.record.num_samplings):
                    vlr.sampling_records[x] =SamplingRecord(pulsebinary)   
                                    
            #if VLR not a scanner or pulse descriptor just read data but do not parse
            #TODO: add additional vlr types                        
            else:        
                vlr.record=	"".join(struct.unpack("c"*vlr.record_length, pulsebinary.read(vlr.record_length))).strip("\x00")
            
            #add vlr to the vlrs dictionary
            self.vlrs[vlr.record_id] = vlr
    
        #close the pls file
        pulsebinary.close()
        
    def get_pulse(self,pulse_number):
        """Given pulse number(a) return the corresponding pulse record(s)
        :param pulse_record: Int or list of pulse number or a pulse record object
        """
        
        #check if pulse number if within range of expected number of pulse
        if pulse_number > self.num_pulses or pulse_number <0:
            print "ERROR: Pulse number outside the range of expected values"
            return

        #initialize pulse object
        pulsebinary =  open(self.filename, 'rb')
        
        record = PulseRecord(pulsebinary,pulse_number,self)
        
        pulsebinary.close()

        return record
                
    def get_waves(self,pulse_record, filename = None):    
        """Give pulse record(s) or pulse number(s) this functions return the corresponding waves
        
        :param pulse_record: Int, pulse number or a pulse record object   
        :param filename: String, pathname to uncompressed waveform file (*.wvs) if non is specified assumes
                         that the waveform used the same pathname as the pulsewaves file but with a ".wvs" extension
        
        """
        
        #if a pulse number is given get the corresponding pulse record
        if type(pulse_record) == int:
            pulse_record = self.get_pulse(pulse_record)
        
        #if input is not a pulse record
        elif type(pulse_record).__name__ == "pulse":
            print "Unrecognized pulse input type, enter pulse number or pulse record"
            return

        wave = Waves(self,pulse_record)
        
        return wave
        
    def print_table(self):
        for key, value in sorted(self.__dict__.iteritems()):
            if type(value) <> dict:
                print "{:<20} {:<15}".format(key, value)
            elif type(value) == dict:
                print "{:<20} {:<15}".format(key, value.keys())
                
                
    def cycle_pulses(self, start, end):
        """Yields a generator that cycles through pulses
           :param start: Int, starting pulse number, default: 0
           :param end: Int, ending pulse number, default: last pulse
           """
        if start >= end:
            print "ERROR: Starting pulse number greater than ending pulse number"
            return
        
        pulsebinary =  open(self.filename, 'rb')
        for pulse_number in range(start,end):
            record = PulseRecord(pulsebinary,pulse_number,self)
        
            yield record

        pulsebinary.close()

    def create_spatial_index(self, overwrite = False):
        """Create and 2D RTree spatial index using last sample coordinates"""
         
        #check if spatial index exists
         
        spatial_index = index.Index(os.path.splitext(self.filename)[0])
         
        print "Generating spatial index.....%s points...this may take a while...." % self.num_pulses
         
        for pulse_record in self.cycle_pulses(0,self.num_pulses):           
            x_last = pulse_record.x_anchor + pulse_record.last_return * pulse_record.dx
            y_last = pulse_record.y_anchor + pulse_record.last_return * pulse_record.dy
            #change to int to save space
            x_last = int(x_last)
            y_last = int(y_last)
            
            #add pulse to index
            spatial_index.insert(pulse_record.pulse_number,(x_last,y_last))
         
            #print status update
            if pulse_record.pulse_number in (np.linspace(0,1,11) * self.num_pulses).astype(int):
                print "%s percent complete." % int(100*pulse_record.pulse_number/float(self.num_pulses))
         
        spatial_index.close() 

    def get_spatial_points(self,x,y,distance):
        """Use spatial index to retrieve pulsewave within a given bounding box
           :param left:
           :param right:
           :param bottom:
           :param top:     
        """
                
        if not os.path.isfile(os.path.splitext(self.filename)[0] + ".idx"):
            print "Spatial index not found!!!"
            return 

        
        spatial_index = index.Index(os.path.splitext(self.filename)[0])

        intersected_pulses = list(spatial_index.intersection((x-distance,y-distance,x+distance,y+distance)))

        spatial_index.close() 


        return intersected_pulses

def openPLS(filename):  
    """Open an uncompressed pulsewaves files (*.pls)
       :param filename: pulsewaves file path
 
    """
    
    return PulseWaves(filename)

          
class PulseRecord(object):
    """Pulse descriptor object"""
 
    def __init__(self,pulsebinary,pulse_number,header):
                
        #jump to the start of the pulse record
        pulsebinary.seek(header.offset_to_pulses+ pulse_number * header.pulse_size)
        
        self.gps_timestamp = header.t_scale * struct.unpack("q", pulsebinary.read(8))[0] + header.t_offset
        self.offset_to_waves = struct.unpack("q", pulsebinary.read(8))[0]
        self.x_anchor =	header.x_scale * struct.unpack("=l", pulsebinary.read(4))[0] + header.x_offset
        self.y_anchor = header.y_scale * struct.unpack("=l", pulsebinary.read(4))[0] + header.y_offset
        self.z_anchor = header.z_scale * struct.unpack("=l", pulsebinary.read(4))[0] + header.z_offset
        self.x_target =	header.x_scale * struct.unpack("=l", pulsebinary.read(4))[0] + header.x_offset
        self.y_target =	header.y_scale * struct.unpack("=l", pulsebinary.read(4))[0] + header.y_offset
        self.z_target =	header.z_scale * struct.unpack("=l", pulsebinary.read(4))[0] + header.z_offset
        self.first_return =  struct.unpack("h", pulsebinary.read(2))[0]
        self.last_return = 	 struct.unpack("h", pulsebinary.read(2))[0]
        self.pulse_number = pulse_number
        bits = []
        for bit in (ord(bit) for bit in pulsebinary.read(2)):
            for i in xrange(8):
                bits.append((bit >> i) & 1)
                                      
        self.pulse_descriptor = 	200000 + int("".join(str(x) for x in bits[:8][::-1]),2)
        self.reserved = bits[8:12]
        self.edge = 	bits[12]
        self.scan_direction = 	bits[13]
        self.facet = 	bits[14:]
        self.intensity = struct.unpack("B", pulsebinary.read(1))[0]
        self.classification = 	struct.unpack("B", pulsebinary.read(1))[0]
        
        #calculate direction vector
        self.dx = (self.x_target - self.x_anchor)/1000
        self.dy = (self.y_target - self.y_anchor)/1000
        self.dz = (self.z_target - self.z_anchor)/1000
        

    def print_table(self):
        for key, value in sorted(self.__dict__.iteritems()):
            if type(value) <> dict:
                print "{:<20} {:<15}".format(key, value)
            elif type(value) == dict:
                print "{:<20} {:<15}".format(key, value.keys())
          
        
class Waves(object):
    
    def __init__(self, header, pulse_record):    
        
        #get sample records corresponding to the waveform
        descriptor = header.vlrs[pulse_record.pulse_descriptor]
        sample_records = descriptor.sampling_records
             
        #read header
        self.filename = os.path.splitext(header.filename)[0]+'.wvs'
        wavebinary =  open(self.filename, 'rb')

        self.file_sig ="".join(struct.unpack("c"*16, wavebinary.read(16))).strip("\x00")
        self.compression = struct.unpack("I", wavebinary.read(4))[0]
        self.reserved = struct.unpack("B"*40, wavebinary.read(40))
        self.segments= {}
        
        #jump to the start of the wave
        wavebinary.seek(pulse_record.offset_to_waves)

        #cycle through each sample
        for key in sample_records.keys():

            sample_record = sample_records[key]           
            duration_anchor = struct.unpack("=L", wavebinary.read(sample_record.bits_anchor/8))[0]   
            num_samples = struct.unpack("=h", wavebinary.read(sample_record.bits_samples/8))[0] 
        
            samples= []
            for sample_num in range(num_samples): 
                          
                sample = struct.unpack("=h", wavebinary.read(sample_record.bits_per_sample/8))[0]
                #calculate 3 dimensional sample coordinates
                x = pulse_record.x_anchor + (duration_anchor + sample_num)  * pulse_record.dx
                y  =pulse_record.y_anchor + (duration_anchor + sample_num)  * pulse_record.dy
                z  =pulse_record.z_anchor + (duration_anchor + sample_num)  * pulse_record.dz
                samples.append([x,y,z,sample])
            
            self.segments[key] = np.array(samples).T
        
        wavebinary.close()
        
        
    def plot(self,save_path = None):
        '''Simple waveform plotting function waveform vs z
        :param segment: Int, segement number to plot
        '''
        
        fig = plt.figure(figsize = (6,4))
        ax_outgoing = fig.add_subplot(121)
        ax_returning = fig.add_subplot(122)        
        ax_outgoing.plot(wave.segments[0][3],wave.segments[0][2],'b')

        ax_outgoing.set_xlabel(r"$\mathrm{\mathsf{Outgoing\/waveform}}$", fontsize= 12)
        ax_outgoing.set_ylabel(r"$\mathrm{\mathsf{Z}}$", fontsize= 12)
        ax_returning.plot(wave.segments[1][3],wave.segments[1][2],'r')
        ax_returning.set_xlabel(r"$\mathrm{\mathsf{Returning\/waveform}}$", fontsize= 12)
        
        #tweak axis settings
        for ax in [ax_outgoing,ax_returning]:
            #hide right and top axes
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            #change axes line widths
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            #change tick font size
            ax.tick_params(labelsize=10)
     
        if save_path:
            plt.savefig(save_path,dpi=600, bbox_inches="tight")
        plt.show()
        plt.close()
     
    def smooth(self, window,polyorder,deriv = 0):
        ''' Use savitzky golay filter to smooth the waveform inplace
        :param window: Int, odd numbered window size
        :param polyorder: Int, smoothing polymomial window
        :param deriv: Int, derivative number
        See scipy.signal.savgolfilter docs for more information
        '''
        #cycle through the segments
        for key,value in self.segments.iteritems():
            self.segments[key][3] = savgol_filter(self.segments[key][3],window,polyorder,deriv)
        
    def print_table(self):
        for key, value in sorted(self.__dict__.iteritems()):
            if type(value) <> dict:
                print "{:<20} {:<15}".format(key, value)
            elif type(value) == dict:
                print "{:<20} {:<15}".format(key, value.keys())
        

class VLR(object):
    """Variable length record (VLR) object"""  
    
    def __init__(self,pulsebinary):    
        self.user_id =	 "".join(struct.unpack("c"*16, pulsebinary.read(16))).strip("\x00").strip("\x00")
        self.record_id =   struct.unpack("I", pulsebinary.read(4))[0]
        self.reserved =	struct.unpack("I", pulsebinary.read(4))[0]
        self.record_length =	struct.unpack("q", pulsebinary.read(8))[0]
        self.desciption = "".join(struct.unpack("c"*64, pulsebinary.read(64))).strip("\x00").strip("\x00")

    def print_table(self):
        for key, value in sorted(self.__dict__.iteritems()):
            if type(value) <> dict:
                print "{:<20} {:<15}".format(key, value)
            elif type(value) == dict:
                print "{:<20} {:<15}".format(key, value.keys())

class Scanner(object):
    
    def __init__(self,pulsebinary):  
        self.size = struct.unpack("I", pulsebinary.read(4))[0]
        self.reserved =	struct.unpack("I", pulsebinary.read(4))[0]
        self.instrument = "".join(struct.unpack("c"*64, pulsebinary.read(64))).strip("\x00").strip("\x00")
        self.serial = "".join(struct.unpack("c"*64, pulsebinary.read(64))).strip("\x00").strip("\x00")
        self.wavelength = struct.unpack("f", pulsebinary.read(4))[0]
        self.out_pulse_width = struct.unpack("f", pulsebinary.read(4))[0]
        self.scan_pattern = struct.unpack("I", pulsebinary.read(4))[0]
        self.num_facets = struct.unpack("I", pulsebinary.read(4))[0]
        self.scan_frequency = struct.unpack("f", pulsebinary.read(4))[0]	 
        self.scan_angle_min = struct.unpack("f", pulsebinary.read(4))[0]	 
        self.scan_angle_max = struct.unpack("f", pulsebinary.read(4))[0]	 
        self.pulse_frequency = struct.unpack("f", pulsebinary.read(4))[0] 
        self.beam_diam = struct.unpack("f", pulsebinary.read(4))[0]
        self.beam_diverge = struct.unpack("f", pulsebinary.read(4))[0]
        self.min_range = struct.unpack("f", pulsebinary.read(4))[0]
        self.max_range = struct.unpack("f", pulsebinary.read(4))[0]
        self.description = "".join(struct.unpack("c"*64, pulsebinary.read(64))).strip("\x00").strip("\x00")

    def print_table(self):
        for key, value in sorted(self.__dict__.iteritems()):
            if type(value) <> dict:
                print "{:<20} {:<15}".format(key, value)
            elif type(value) == dict:
                print "{:<20} {:<15}".format(key, value.keys())
                
class SamplingRecord(object):
    
    
    def __init__(self,pulsebinary):     
        
        self.size = struct.unpack("I", pulsebinary.read(4))[0]
        self.reserved =	struct.unpack("I", pulsebinary.read(4))[0]
    
        sample_type = {1:"outgoing",2: "returning"}
        self.type = sample_type[struct.unpack("B", pulsebinary.read(1))[0]]
        
        self.channel = struct.unpack("B", pulsebinary.read(1))[0]
        self.unused = struct.unpack("B", pulsebinary.read(1))[0]
        self.bits_anchor = struct.unpack("B", pulsebinary.read(1))[0]
        self.scale_anchor = struct.unpack("f", pulsebinary.read(4))[0]
        self.offset_anchor = struct.unpack("f", pulsebinary.read(4))[0]
        self.bits_segments = struct.unpack("B", pulsebinary.read(1))[0]
        self.bits_samples = struct.unpack("B", pulsebinary.read(1))[0]
        self.num_segments = struct.unpack("H", pulsebinary.read(2))[0]
        self.num_samples =  struct.unpack("I", pulsebinary.read(4))[0]
        self.bits_per_sample = struct.unpack("H", pulsebinary.read(2))[0]
        self.lut_index = struct.unpack("H", pulsebinary.read(2))[0]
        self.samples_units = struct.unpack("f", pulsebinary.read(4))[0]
        self.compression =  struct.unpack("I", pulsebinary.read(4))[0]
        self.description = "".join(struct.unpack("c"*64, pulsebinary.read(64))).strip("\x00").strip("\x00")   

    def print_table(self):
        for key, value in sorted(self.__dict__.iteritems()):
            if type(value) <> dict:
                print "{:<20} {:<15}".format(key, value)
            elif type(value) == dict:
                print "{:<20} {:<15}".format(key, value.keys())

class PulseDecriptor(object):

    def __init__(self,pulsebinary):  
        self.size = struct.unpack("I", pulsebinary.read(4))[0]
        self.reserved =	struct.unpack("I", pulsebinary.read(4))[0]
        self.optical_center = struct.unpack("=l", pulsebinary.read(4))[0]
        self.num_extra_wave_bytes = struct.unpack("H", pulsebinary.read(2))[0]
        self.num_samplings = struct.unpack("H", pulsebinary.read(2))[0]
        self.sample_units = struct.unpack("f", pulsebinary.read(4))[0]
        self.compression = struct.unpack("I", pulsebinary.read(4))[0]
        self.scanner_index = struct.unpack("I", pulsebinary.read(4))[0]
        self.description = "".join(struct.unpack("c"*64, pulsebinary.read(64))).strip("\x00").strip("\x00")


    def print_table(self):
        for key, value in sorted(self.__dict__.iteritems()):
            if type(value) <> dict:
                print "{:<20} {:<15}".format(key, value)
            elif type(value) == dict:
                print "{:<20} {:<15}".format(key, value.keys())
        



    
    
    