import socket
import logging
import time
import struct
import io

import numpy as num

from pyrocko import trace

logger = logging.getLogger('UDPStream')
TRACE_BUFFER = 256*3*8


class UDPStream(object):

    def __init__(self, mcast_grp, port):
        self.mcast_grp = mcast_grp
        self.port = int(port)

        self.socket = None
        self.running = False

    def acquisition_start(self):
        assert not self.running
        logger.info('Starting listening to UDP Multicast group %s:%s'
                    % (self.mcast_grp, self.port))
        self.sock = socket.socket(
            family=socket.AF_INET,
            type=socket.SOCK_DGRAM,
            proto=socket.IPPROTO_UDP)

        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind((self.mcast_grp, self.port))
        mreq = struct.pack("4sl", socket.inet_aton(self.mcast_grp),
                           socket.INADDR_ANY)
        self.sock.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)

        self.running = True

    def acquisition_stop(self):
        self.sock.close()

    def process(self):
        print('Receiving data')
        buff = io.BytesIO()
        while buff.tell() < 1024*3:
            buff.write(self.sock.recv(512))

        buff.seek(0)
        data = num.frombuffer(buff.read(), dtype=num.float64)
        # data = data.reshape(data.size//3, 3)

        for icha, cha in enumerate(('SHX')):
            tr = trace.Trace(
                network='UD',
                station='S1',
                location='NT',
                channel=cha,
                tmin=time.time(),
                deltat=1./100,
                ydata=data)
            self.got_trace(tr)

        time.sleep(1.)

        return True

    def got_trace(self, tr):
        logger.info('Got trace from UDP socket: %s' % tr)
