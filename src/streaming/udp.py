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
            self.queue.put_nowait(data)

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
        hdr = comp, nsamples, sampling, t0 = struct.unpack('chhd', data[:16])
        data = data[16:]
        return hdr, data

    def process(self):
        traces = []
        if self.udpqueue.qsize() < 40:
            return True

        for ip in range(self.udpqueue.qsize()):
            traces.append(self.extract_trace(self.udpqueue.get()))

        channels = set([tr[0][0] for tr in traces])
        for cha in channels:
            trs = [tr for tr in traces if tr[0][0] == cha]
            comp, _, sampling, t0 = trs[0][0]

            data = b''.join([tr[1] for tr in trs])
            data = num.frombuffer(data)

            tr = trace.Trace(
                network='UD',
                station='S1',
                location='NT',
                channel='SH' + comp.decode().upper(),
                tmin=t0,
                deltat=1./sampling,
                ydata=data)
            self.got_trace(tr)

        return True

    def got_trace(self, tr):
        log.info('Got trace from UDP socket: %s' % tr)
