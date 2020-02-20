#!/usr/bin/env python3
# -*- mode: python -*-

# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You can redistribute it and/or modify it under
# the terms of the Do What The Fuck You Want To Public License, Version 2, as
# published by Sam Hocevar. See http://sam.zoy.org/wtfpl/COPYING for more
# details.

# Some useful resources:
# - http://atomicparsley.sourceforge.net/mpeg-4files.html
# - http://developer.apple.com/library/mac/#documentation/QuickTime/QTFF/QTFFChap2/qtff2.html
# - http://www.sno.phy.queensu.ca/~phil/exiftool/TagNames/QuickTime.html

import datetime
import argparse
import os.path
import struct
import sys
import time

NAMES = { # not used for anything, but just documents a short blurb
          # about what these things mean
    "vmhd": "video information media header",
    "mvhd": 'movie header',
    "tkhd": 'track header',
    "mdhd": 'media header', # The media header atom specifies the characteristics of a media, including time scale and duration
    "smhd": 'sound media information header', 
    "hdlr": 'handler reference', # The handler reference atom specifies the media handler component that is to be used to interpret the media’s data

    "stsd": "sample description", # The sample description atom contains a table of sample descriptions
    "stts": "time-to-sample", # Time-to-sample atoms store duration information for a media’s samples, providing a mapping from a time in a media to the corresponding data sample
    "stsc": "sample-to-chunk", # The sample-to-chunk atom contains a table that maps samples to chunks
    "stco": 'chunk offset', # Chunk offset atoms identify the location of each chunk of data
    "stsz": 'sample size', # You use sample size atoms to specify the size of each sample
    "ctts": 'composition offset', # The composition offset atom contains a sample-by-sample mapping of the decode-to-presentation time
    "stss": "sync sample", # The sync sample atom identifies the key frames
}

parser = argparse.ArgumentParser(description='Parse and modify quicktime atoms.')
parser.add_argument('mov_file', metavar='file.mov')
parser.add_argument('--mdat_length', metavar='N', type=int)
parser.add_argument('--number_of_samples', metavar='M', type=int)
parser.add_argument('--chunk_table_file', type=argparse.FileType('rb'))
parser.add_argument('--sample_size_table', type=argparse.FileType('br'))
parser.add_argument('--out_moov_file', type=argparse.FileType('bw+'))
parser.add_argument('--new_moov_file', type=argparse.FileType('br'))


CONTAINER_ATOMS = ["moov", "trak", "mdia", "minf","dinf","stbl"]
_IGNORE_ATOMS = ["iods"] # couldn't find documentation for this
_ATOMS = {
    "pnot": (12, "I2x4s2x",
             ("Modification time", "Atom type"),
             (0,)),
    "vmhd": (12, "4xH6x",
             ("graphics mode",),
             () ),
    "mvhd": (100, "4x5IH10x36x7I",
             ("Creation time", "Modification time",
              "Time Scale",
              'Duration',
              'Preferred rate',
              'Preferred volume',
              'preview time',
              'preview duration',
              'poster time',
              'selection time',
              'selection duration',
              'current time',
              'next track id'
              ),
             (4, 8)),
    "tkhd": (84, "4x2I72x",
             ("Creation time", "Modification time"),
             (4, 8)),
    "mdhd": (24, "B3x4I2H", #3x is "flags"
             ("Version", "Creation time", "Modification time","Time Scale","Duration","Language","Quality"),
             (4, 8)), # positions where dates are so we can modify them
    "smhd": (8, "4xH2x",
             ("balance",),
             ())
}

_VARIABLE_LEN_ATOMS = {
    "hdlr": (4 + 5*4, "4x5I",
             ("Component type",
              'component subtype',
              'component manufacturer',
              'component flags',
              'component flags mask'),
             (),
             "component name"
             ),
    "stsd": (8, "4xI",
             ("number of entries",),
             (),
             "sample description table"),
    "stts": (8,"4xI",
             ("number of entries",),
             (),
             "time-to-sample table"),
    "stsc": (8,"4xI",
             ("number of entries",),
             (),
             "sample-to-chunk table"),
    "stco": (8,"4xI",
             ("number of entries",),
             (),
             "chunk offset table"),
    "stsz": (12,"4xII",
             ("sample size","number of entries",),
             (),
             "sample size table"),
    "ctts": (12,"4xII",
             ("entry count",),
             (),
             "composition offset table"),
    "stss": (12,"4xII",
             ("number of entries",),
             (),
             "sync sample table")
    

}

_VARIABLE_CHAINED_ATOMS = {
    "dref": (8, "4xI",
             ("number of entries",),
             (),
             "data references"
             )
}

_DATES = ("Creation time", "Modification time")

class Mov(object):
    def __init__(self, fn, number_of_samples, chunk_table_file, sample_size_table, moov_out_file):
        self._fn = fn
        self._offsets = []
        self._mdat_offset = -1
        self._moov_offset = -1
        self._trak_offset = -1
        self._number_of_samples = number_of_samples
        self._chunk_table_file = chunk_table_file
        self._sample_size_table = sample_size_table
        self._moov_out_file = moov_out_file
        if self._moov_out_file:
          self._write_moov = True
        else:
          self._write_moov = False

    def parse(self):
        fsize = os.path.getsize(self._fn)
        print("File: {} ({} bytes, {} MB)".format(self._fn, fsize, fsize / (1024.**2)))
        with open(self._fn, "rb") as self._f:
            self._parse(fsize)


    def _f_seek(self,l, stop_write_out_override=False):
        print('seeking '+str(l))
        if self._write_moov and not stop_write_out_override:
          dat = self._f.read(l)
          self._moov_out_file.write(dat)
        else:
          self._f.seek(l, 1)

    def _f_read(self,l, stop_write_out_override=False):
        print('reading '+str(l))
        dat = self._f.read(l)
        if self._write_moov and not stop_write_out_override:
          self._moov_out_file.write(dat)
        return dat


    def _parse(self, length, depth=0):
        prefix = "  "*depth + "- "

        if self._moov_out_file:
            if depth > 0:
              current_container_offset = self._moov_out_file.tell() - 8
            else:
              current_container_offset = 0
        n = 0
        while n < length:
            data = self._f_read(8)
            #print(len(data), data)
            al, an = struct.unpack(">I4s", data)
            an = an.decode()
            print("{}Atom: {} ({} bytes)".format(prefix, an, al))

            if an == "moov":
              self._moov_offset = self._f.tell() - 8

            if an == "trak" and self._moov_out_file:
              self._trak_offset = self._moov_out_file.tell() - 8
              if not self._write_moov:
                self._write_moov = True
                # re-read header so its writtent ofile
                self._f.seek(-8, 1)
                self._f_seek(8)

            if an == "mvhd" and self._write_moov:
                self._write_mvhd(al-8, depth)
            elif an == "mdhd" and self._write_moov:
                self._write_mdhd(al-8, depth)
            elif an in _ATOMS:
                self._parse_atom(an, al-8, depth)
            elif an == "udta":
                self._parse_udta(al-8, depth)
            elif an == "stts":
              if self._write_moov:
                self._write_stts(al-8, depth)
              else:
                self._parse_stts(al-8, depth)
            elif an == "stsc":
              if self._write_moov and self._last_component_subtype == b'vide':
                self._write_stsc(al-8, depth)
              else:
                self._parse_stsc(al-8, depth)
            elif an == "stsz":
                # read and rewind header
                dat = self._f_read(8, stop_write_out_override=True)
                sample_size = struct.unpack(">4xI", dat)[0]
                self._f.seek(-8, 1)
                if self._sample_size_table and self._write_moov and sample_size == 0:
                  self._write_out_stsz(al-8, depth)
                else:
                  self._parse_stsz(al-8, depth)
            elif an == "hdlr":
                self._parse_hdlr(al-8, depth)
            elif an == "co64":
                if self._chunk_table_file and self._write_moov and self._last_component_subtype == b'vide':
                  self._write_out_co64(al-8, depth)
                else:
                  self._parse_co64(al-8, depth)
            elif an == "ftyp":
                self._read_ftyp(al-8, depth)
            elif an in CONTAINER_ATOMS:
                self._parse(al-8, depth+1)

            elif an in _VARIABLE_LEN_ATOMS:
                self._parse_atom(an, al-8, depth, variable=True)
            elif an in _VARIABLE_CHAINED_ATOMS:
                self._parse_atom(an, al-8, depth, chained=True)
            elif an in _IGNORE_ATOMS:
                self._f_read(al-8)
            else:
                print('unhandled thingie',al,an)
                loc = self._f.tell() - 8
                if al == 1:
                    # 64 bit!
                    print("64 bit header!")
                    al = struct.unpack(">Q", self._f_read(8))[0]
                    long_header = True
                    self._f_seek(al-16)
                elif al == 0 and depth == 0:
                  print("atom goes to EOF")
                  al = length - n
                  self._f_seek(al-8)
                  long_header = False
                else:
                    self._f_seek(al-8)
                    long_header = False

                if an == "mdat":
                  self._mdat_offset = loc
                  self._mdat_64_loc = long_header
                  current_byte= self._f.tell()
                  if long_header:
                    self._f.seek(loc+16, 0)
                  else:
                    self._f.seek(loc+8, 0)
                  m = 0
                  while m < 2:
                    sub_length, subtype = struct.unpack(">I4s", self._f_read(8))
                    print("  - {}{}: {} - {}".format(prefix, m, sub_length, subtype))
                    self._f.seek(sub_length - 8, 1)
                    m = m + 1
                  raw = self._f_read(12)
                  print("  - {}raw: {}".format(prefix, raw))


                  self._f.seek(current_byte, 0)

            n += al
        if self._write_moov:
            current_pos = self._moov_out_file.tell()
            length = current_pos - current_container_offset
            print("{}rewriting container length to {} bytes, offset: {}, cp:{}".format(prefix, length, current_container_offset, current_pos))
            self._moov_out_file.seek(current_container_offset, 0)
            self._moov_out_file.write(struct.pack(">I", length))
            self._moov_out_file.seek(current_pos, 0)


    def _parse_atom(self, atom, length, depth, variable=False, chained=False):
        if variable:
            spec = _VARIABLE_LEN_ATOMS[atom]
        elif chained:
            spec = _VARIABLE_CHAINED_ATOMS[atom]
        else:
            spec = _ATOMS[atom]
            assert length == spec[0]
            
        pos = self._f.tell()
        prefix = "  "*depth + "  | "
        data = self._f_read(length)
        if variable:
            v = struct.unpack(">"+spec[1], data[:spec[0]])
        elif chained:
            v = struct.unpack(">"+spec[1], data[:spec[0]])
        else:
            v = struct.unpack(">"+spec[1], data)
        k = spec[2]
        for i in range(len(k)):
            vv = v[i]
            if type(vv) == bytes:
                vv = vv.decode()
            elif k[i] in _DATES:
                vv = self._macdate2date(vv)
            print("{}{}: {}".format(prefix, k[i], vv))

        if variable or chained:
            lim = 100
            realdata = data[spec[0]:]
            if len(realdata) > lim:
                print("{}{}: {}{}{}{}".format(prefix, spec[4], realdata[:lim], '...', len(realdata)-lim,' more bytes'))
            else:
                print("{}{}: {}".format(prefix, spec[4], realdata))

        for offset in spec[3]:
            self._offsets.append(pos + offset)

    def _read_ftyp(self, length, depth):
        prefix = "  "*depth + "  | "
        data = self._f_read(8)
        brand, version = struct.unpack(">4sI", data)
        brand = brand.decode("latin1")
        print("{}Brand: {}, version: {}".format(prefix, brand, version))
        self._f_read(length-8)

    def _parse_udta(self, length, depth):
        prefix = "  "*depth + "  | "
        n = 0
        while n < length:
            atom_size, data_type = struct.unpack(">I4s", self._f_read(8))
            data_type = data_type.decode("latin1")
            if atom_size == 0:
              print("skip 0 length")
              break
            raw = self._f_read(atom_size-8)
            if data_type[0] == "©":
                print("{}{}: {}".format(prefix, data_type, raw[3:].decode()))
            else:
                print("{}{} ({} bytes) {}".format(prefix, data_type, atom_size-8, raw.decode()))
            n += atom_size

    def _write_mvhd(self, length, depth):
        prefix = "  "*depth + "  | "
        # the length was already written out, so rewind
        print("{}rewinding moov out 8 bytes to write out new header".format(prefix))
        self._moov_out_file.seek(-8, 2)
        
        # header
        # 4 - atom size
        # 4 - atom type - mvhd
        # 16 - more header
        # 4 - duration
        # 76 - more header
        # 4 - next track id

        header = struct.pack(">I4s", 108, b'mvhd')
        # copy Version/flags over from original
        header = header + self._f_read(16, stop_write_out_override=True)
        length = length - 16
        header = header + struct.pack(">I", self._number_of_samples*1001)
        self._f_seek(4, stop_write_out_override=True)
        header = header + self._f_read(76, stop_write_out_override=True)
        # we delete the audio track so reduce next track id
        header = header + struct.pack(">I", 3)
        self._f_seek(4, stop_write_out_override=True)
        self._moov_out_file.write(header)

        data = header[8:]
        spec = _ATOMS["mvhd"]
        v = struct.unpack(">"+spec[1], data)
        k = spec[2]
        for i in range(len(k)):
            vv = v[i]
            if type(vv) == bytes:
                vv = vv.decode()
            elif k[i] in _DATES:
                vv = self._macdate2date(vv)
            print("{}{}: {}".format(prefix, k[i], vv))


    def _write_mdhd(self, length, depth):
        prefix = "  "*depth + "  | "
        # the length was already written out, so rewind
        print("{}rewinding moov out 8 bytes to write out new header".format(prefix))
        self._moov_out_file.seek(-8, 2)
        
        # header
        # 4 - atom size
        # 4 - atom type - mdhd
        # 16 - more header
        # 4 - duration
        # 4 - more header

        header = struct.pack(">I4s", 32, b'mdhd')
        # copy Version/flags over from original
        header = header + self._f_read(16, stop_write_out_override=True)
        length = length - 16
        duration = struct.unpack(">I", self._f_read(4, stop_write_out_override=True))[0]
        multiple = int(duration/7602) # special number, number of samples in template moov file

        print("{}multiples {} - og duratin {}".format(prefix, multiple, duration))
        header = header + struct.pack(">I", self._number_of_samples*multiple)
        header = header + self._f_read(4, stop_write_out_override=True)
        self._moov_out_file.write(header)

        data = header[8:]
        spec = _ATOMS["mdhd"]
        v = struct.unpack(">"+spec[1], data)
        k = spec[2]
        for i in range(len(k)):
            vv = v[i]
            if type(vv) == bytes:
                vv = vv.decode()
            elif k[i] in _DATES:
                vv = self._macdate2date(vv)
            print("{}{}: {}".format(prefix, k[i], vv))

    def _parse_stts(self, length, depth):
        prefix = "  "*depth + "  | "
        dat = self._f_read(8)
        num_entries = struct.unpack(">I", dat[4:])[0]
        print("{}number of entries : {}".format(prefix, num_entries))
        print("{}time-to-sample table:".format(prefix))
        n = 0
        while n < num_entries and n < 10:
            sample_count, sample_duration = struct.unpack(">II", self._f_read(8))
            print("{}{}: {}".format(prefix, sample_count, sample_duration))
            n = n + 1

    def _write_stts(self, length, depth):
        prefix = "  "*depth + " writing | "

        stts_raw = struct.pack(">II", self._number_of_samples, 1001)
        # add 16 for the header length, 8 for the above
        stts_table_size = 16 + 8

        # the length was already written out, so rewind
        print("{}rewinding moov out 8 bytes to write out new header".format(prefix))
        self._moov_out_file.seek(-8, 2)
        
        # header
        # 4 - atom size
        # 4 - atom type - stts
        # 1 - Version
        # 3 - Flags
        # 4 - number entries

        header = struct.pack(">I4s", stts_table_size, b'stts')
        # copy Version/flags over from original
        header = header + self._f_read(4, stop_write_out_override=True)
        length = length - 4
        header = header + struct.pack(">I", 1)
        self._moov_out_file.write(header)
        self._moov_out_file.write(stts_raw)
        
        al, an = struct.unpack(">I4s", header[:8])
        an = an.decode()
        print("{}Atom: {} ({} bytes)".format(prefix, an, al))

        dat = header[8:]
        num_entries = struct.unpack(">I", dat[4:])[0]
        print("{}number of entries : {}".format(prefix, num_entries))
        print("{}time-to-sample table:".format(prefix))
        n = 0
        while n < num_entries and n < 10:
            sample_count, sample_duration = struct.unpack(">II", stts_raw[n*8:n*8+8])
            print("{}{}: {}".format(prefix, sample_count, sample_duration))
            n = n + 1
        self._f_seek(length, stop_write_out_override=True)

    def _parse_stsc(self, length, depth):
        prefix = "  "*depth + "  | "
        dat = self._f_read(8)
        num_entries = struct.unpack(">I", dat[4:])[0]
        print("{}number of entries : {}".format(prefix, num_entries))
        print("{}sample-to-chunk table:".format(prefix))
        n = 0
        while n < num_entries and n < 10:
            first_chunk, samples_per_chunk, sample_description_id = struct.unpack(">III", self._f_read(12))
            print("{}{}: {} - {}".format(prefix, first_chunk, samples_per_chunk, sample_description_id))
            n = n + 1


    def _write_stsc(self, length, depth):
        prefix = "  "*depth + " writing | "

        stsc_raw = struct.pack(">III", 1, 1, 1)
        stsc_raw = stsc_raw + struct.pack(">III", self._number_of_samples, 1, 1)
        # add 16 for the header length plus 12 for each entry
        stsc_table_size = 16 + 12*2

        # the length was already written out, so rewind
        print("{}rewinding moov out 8 bytes to write out new header".format(prefix))
        self._moov_out_file.seek(-8, 2)
        
        # header
        # 4 - atom size
        # 4 - atom type - stsz
        # 1 - Version
        # 3 - Flags
        # 4 - number entries

        header = struct.pack(">I4s", stsc_table_size, b'stsc')
        # copy Version/flags size over from original
        header = header + self._f_read(4, stop_write_out_override=True)
        length = length - 4
        header = header + struct.pack(">I", 2)
        self._moov_out_file.write(header)
        self._moov_out_file.write(stsc_raw)
                
        al, an = struct.unpack(">I4s", header[:8])
        an = an.decode()
        print("{}Atom: {} ({} bytes)".format(prefix, an, al))
        
        dat = header[8:]
        num_entries = struct.unpack(">I", dat[4:])[0]
        print("{}number of entries : {}".format(prefix, num_entries))
        n = 0
        while n < num_entries and n < 10:
            first_chunk, samples_per_chunk, sample_description_id = struct.unpack(">III", stsc_raw[n*12:n*12+12])
            print("{}{}: {} - {}".format(prefix, first_chunk, samples_per_chunk, sample_description_id))
            n = n + 1
        self._f_seek(length, stop_write_out_override=True)

    def _parse_stsz(self, length, depth):
        prefix = "  "*depth + "  | "
        dat = self._f_read(12)
        length = length - 12
        sample_size, num_entries= struct.unpack(">II", dat[4:])

        print("{}sample size : {}".format(prefix, sample_size))
        print("{}number of entries : {}".format(prefix, num_entries))

        if sample_size != 0:
          return
        print("{}sample size table:".format(prefix))
        n = 0
        while n < num_entries and n < 10:
            ind_sample_size = struct.unpack(">I", self._f_read(4))[0]
            length = length - 4
            print("{}{} - {}".format(prefix, n, ind_sample_size))
            n = n + 1
        self._f_seek(length)


    def _write_out_stsz(self, length, depth):
        prefix = "  "*depth + " writing | "

        sample_size_raw = self._sample_size_table.read()
        # add 20 for the header length
        stsz_table_size = 20 + self._sample_size_table.tell()

        # the length was already written out, so rewind
        print("{}rewinding moov out 8 bytes to write out new header".format(prefix))
        self._moov_out_file.seek(-8, 2)
        
        # header
        # 4 - atom size
        # 4 - atom type - stsz
        # 1 - Version
        # 3 - Flags
        # 4 - Sample size, 0 in this case since we have a table
        # 4 - number entries

        header = struct.pack(">I4s", stsz_table_size, b'stsz')
        # copy Version/flags/sample size over from original
        header = header + self._f_read(8, stop_write_out_override=True)
        length = length - 8
        header = header + struct.pack(">I", self._number_of_samples)
        self._moov_out_file.write(header)
        self._moov_out_file.write(sample_size_raw)

        al, an = struct.unpack(">I4s", header[:8])
        an = an.decode()
        print("{}Atom: {} ({} bytes)".format(prefix, an, al))
        
        dat = header[8:]
        sample_size, num_entries= struct.unpack(">II", dat[4:])

        print("{}sample size : {}".format(prefix, sample_size))
        print("{}number of entries : {}".format(prefix, num_entries))

        print("{}sample size table:".format(prefix))
        n = 0
        while n < num_entries:
            ind_sample_size = struct.unpack(">I", sample_size_raw[4*n:4*n+4])[0]
            print("{}{} - {}".format(prefix, n, ind_sample_size))
            n = n + 1
        self._f_seek(length, stop_write_out_override=True)


    def _parse_hdlr(self, length, depth):
        prefix = "  "*depth + "  | "
        dat = self._f_read(24)
        length = length - 24
        component_type, component_subtype, manufacturer, flags, flag_masks = struct.unpack(">4x4s4s3I", dat)
        print("{}Component type: {}".format(prefix, component_type))
        print("{}component subtype: {}".format(prefix, component_subtype))
        print("{}component manufacturer: {}".format(prefix, manufacturer))
        print("{}component flags: {}".format(prefix, flags))
        print("{}component flags mask: {}".format(prefix, flag_masks))

        if component_type == b'mhlr':
          self._last_component_subtype = component_subtype

        if component_subtype == b'soun' and self._write_moov:

          print("{}component subtype == soun, erasing back to {}".format(prefix, self._trak_offset))
          self._write_moov = False
          self._moov_out_file.seek(self._trak_offset, 0)
          self._moov_out_file.truncate()
          self._trak_offset =  -1


        dat = self._f_read(1)
        length = length - 1
        component_name_length = ord(dat)
        print("{}component name length: {}".format(prefix, component_name_length))

        dat = self._f_read(component_name_length)
        length = length - component_name_length
        component_name = struct.unpack(">{}s".format(component_name_length), dat)[0]
        print("{}component name: {}".format(prefix, component_name))

        self._f_seek(length)

    def _parse_co64(self, length, depth):
        prefix = "  "*depth + "  | "
        dat = self._f_read(8)
        length = length - 8
        num_entries = struct.unpack(">I", dat[4:])[0]
        print("{}number of entries : {}".format(prefix, num_entries))
        print("{}chunk offset table:".format(prefix))
        n = 0
        while n < num_entries and n < 10:
            chunk_offset = struct.unpack(">Q", self._f_read(8))[0]
            length = length - 8
            print("{}{}: {}".format(prefix, n, chunk_offset))

            # current_byte= self._f.tell()
            # self._f.seek(chunk_offset, 0)
            # m = 0
            # while m < 2:
            #   sub_length, subtype = struct.unpack(">I4s", self._f_read(8))
            #   print("  - {}{}: {} - {}".format(prefix, m, sub_length, subtype))
            #   self._f.seek(sub_length - 8, 1)
            #   m = m + 1
            # #raw = self._f_read(687)
            # #print("  - {}raw: {}".format(prefix, raw))


            # self._f.seek(current_byte, 0)

            n = n + 1
        self._f_seek(length)

    def _write_out_co64(self, length, depth):
        prefix = "  "*depth + " writing | "

        chunk_table_raw = self._chunk_table_file.read()
        # add 16 for the header length
        chunk_table_size = 16 + self._chunk_table_file.tell()

        # the length was already written out, so rewind
        print("{}rewinding moov out 8 bytes to write out new header".format(prefix))
        self._moov_out_file.seek(-8, 2)
        
        # header
        # 4 - atom size
        # 4 - atom type - co64
        # 1 - Version
        # 3 - Flags
        # 4 - number entries

        header = struct.pack(">I4s", chunk_table_size, b'co64')
        # copy Version/flags over from original
        header = header + self._f_read(4, stop_write_out_override=True)
        length = length - 4
        header = header + struct.pack(">I", self._number_of_samples)
        self._moov_out_file.write(header)
        self._moov_out_file.write(chunk_table_raw)
        
        al, an = struct.unpack(">I4s", header[:8])
        an = an.decode()
        print("{}Atom: {} ({} bytes)".format(prefix, an, al))

        dat = header[8:]
        num_entries = struct.unpack(">I", dat[4:])[0]
        print("{}number of entries : {}".format(prefix, num_entries))
        print("{}chunk offset table:".format(prefix))
        n = 0
        while n < num_entries and n < 10:
            chunk_offset = struct.unpack(">Q", chunk_table_raw[n*8:n*8+8])[0]
            print("{}{}: {}".format(prefix, n, chunk_offset))

            # current_byte= self._f.tell()
            # self._f.seek(chunk_offset, 0)
            # m = 0
            # while m < 2:
            #   sub_length, subtype = struct.unpack(">I4s", self._f_read(8))
            #   print("  - {}{}: {} - {}".format(prefix, m, sub_length, subtype))
            #   self._f.seek(sub_length - 8, 1)
            #   m = m + 1
            # #raw = self._f_read(687)
            # #print("  - {}raw: {}".format(prefix, raw))


            # self._f.seek(current_byte, 0)

            n = n + 1

        self._f_seek(length, stop_write_out_override=True)


    def _macdate2date(self, md):
        d = datetime.datetime(1904, 1, 1) + datetime.timedelta(seconds=md)
        return "{} ({})".format(d, md)

    def _date2macdate(self, d):
        td = datetime.datetime(1970, 1, 1) - datetime.datetime(1904, 1, 1)
        dd = d + td
        sec = time.mktime(dd.timetuple()) - time.timezone
        return int(sec)

    def set_mdat_length(self, l):
        print("New mdat: {}".format(l))
        if self._mdat_offset == -1:
          print("mdat not found in file")
          return
        with open(self._fn, "r+b") as f:
            print("Writing new mdat at {} position".format(self._mdat_offset))
            if self._mdat_64_loc: # the length is stored in the 64 bit format
              f.seek(self._mdat_offset + 8)
              data = struct.pack(">Q", l)
              f.write(data)
            else: # the length is stored in the 32 bit format
              f.seek(self._mdat_offset)
              data = struct.pack(">I", l)
              f.write(data)

            f.flush()
        print("Done! :)")

    def set_moov_data(self, _new_moov_file):
        print("New moov")
        if self._mdat_offset == -1:
          print("mdat not found in file")
          return
        moov_raw = _new_moov_file.read()
        with open(self._fn, "r+b") as f:
            print("Writing new moov at {} position".format(self._moov_offset))

            f.seek(self._moov_offset)
            f.write(moov_raw)
            f.truncate()

            f.flush()
        print("Done! :)")


    def set_date(self, d):
        md = self._date2macdate(d)
        print("New date: {} ({})".format(d, md))
        with open(self._fn, "r+b") as f:
            print("Writing new date at {} positions...".format(len(self._offsets)))
            for offset in self._offsets:
                f.seek(offset)
                data = struct.pack(">I", md)
                f.write(data)
            f.flush()
            print("Touching file...")
            ts = time.mktime(d.timetuple())
            os.utime(self._fn, (ts, ts))
        print("Done! :)")

if __name__ == "__main__":
    usage = "Usage: %prog [options] file.mov [\"YYYY-MM-DD hh:mm:ss\"]"
    args = parser.parse_args()

    m = Mov(args.mov_file, args.number_of_samples, args.chunk_table_file, args.sample_size_table, args.out_moov_file)
    m.parse()
    if args.out_moov_file:
      args.out_moov_file.flush()
      args.out_moov_file.close()

    if args.mdat_length:
        m.set_mdat_length(args.mdat_length)


    if args.new_moov_file:
        m.set_moov_data(args.new_moov_file)
