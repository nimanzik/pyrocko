import socket
import logging
import time
import struct
import threading
import queue

import numpy as num

from pyrocko import trace

log = logging.getLogger('UDPStream')
TRACE_BUFFER = 256*3*8


class UDPCollector(threading.Thread):

    def __init__(self, queue, mcast_grp, port):
        threading.Thread.__init__(self)
        self.mcast_grp = mcast_grp
        self.port = port
        self.queue = queue
        self.socket = self.get_socket(mcast_grp, int(port))

    def run(self):
        print('Start listening to UDP Multicast group %s:%s'
              % (self.mcast_grp, self.port))
        while True:
            data = self.socket.recv(1024)
            self.queue.put(data)

    @staticmethod
    def get_socket(mcast_grp, port):
        sock = socket.socket(
            family=socket.AF_INET,
            type=socket.SOCK_DGRAM,
            proto=socket.IPPROTO_UDP)

        sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        sock.bind((mcast_grp, port))
        mreq = struct.pack("4sl", socket.inet_aton(mcast_grp),
                           socket.INADDR_ANY)
        sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)
        return sock


class UDPStream(object):

    def __init__(self, mcast_grp, port):
        self.mcast_grp = mcast_grp
        self.port = int(port)

        self.running = False
        self.udpqueue = queue.Queue()
        self.collector = UDPCollector(self.udpqueue, mcast_grp, port)

    def acquisition_start(self):
        self.collector.start()

    def acquisition_stop(self):
        self.collector.terminate()

    @staticmethod
    def extract_trace(data):
        header = data[:28]
        nsamples, sampling_rate, t0 = struct.unpack('hhd', header[12:])
        header = {
            'nsl': header[:12],
            'network': header[0:3].decode().strip(),
            'station': header[3:6].decode().strip(),
            'location': header[6:9].decode().strip(),
            'channel': header[9:12].decode().strip(),
            'deltat': 1. / sampling_rate,
            'tmin': t0
        }
        return header, data[28:]

    def process(self):
        traces = []
        if self.udpqueue.qsize() < 20:
            return True

        for ip in range(self.udpqueue.qsize()):
            traces.append(self.extract_trace(self.udpqueue.get()))

        codes = set([tr[0]['nsl'] for tr in traces])
        for code in codes:
            trs = [tr for tr in traces if tr[0]['nsl'] == code]
            data = b''.join([tr[1] for tr in trs])
            data = num.frombuffer(data)

            header = trs[0][0].copy()
            del header['nsl']

            tr = trace.Trace(
                ydata=data,
                **header)
            self.got_trace(tr)

        return True

    def got_trace(self, tr):
        log.info('Got trace from UDP socket: %s' % tr)
