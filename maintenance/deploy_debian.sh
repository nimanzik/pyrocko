#!/bin/bash

GITFORK=pyrocko
GITTARGET=master
# Process command line arguments
while getopts f:t: opt
do
   case "$opt" in
      f) GITFORK=$OPTARG;;
      t) GITTARGET=$OPTARG;;
   esac
done


DEBVERSION=1
DATE=`date +"%a, %d %b %Y %H:%M:%S %z"`

# Setting PATH to correct python distribution, avoid to use virtualenv
export PATH=/usr/bin:/usr/sbin:/bin:/sbin
CODENAME=`lsb_release -cs`
# the lsb-release package in raspbian wheezy
# (https://www.raspberrypi.org/downloads) does not report the codename
# correctly, so fix this
if [ "$CODENAME" == "n/a" ] && [ `arch` == "armv6l" ]; then CODENAME=wheezy; fi
BUILDDIR=/tmp/python-pyrocko_build
PACKAGEDIR=$BUILDDIR/packages
GITDIR=$BUILDDIR/git

# deactivate, else each time all packages are removed
rm -rf $BUILDDIR
mkdir -p $PACKAGEDIR
git clone git@git.pyrocko.org:pyrocko/pyrocko.git $GITDIR
cd $GITDIR
if [ "$GITFORK" != "pyrocko" ]
then
    git remote add upstream git://git.pyrocko.org:pyrocko/pyrocko.git
    git fetch upstream
fi

# Build package
echo "#### Working on $GITTARGET"
cd $GITDIR
git clean -fxd
git fetch --all
git checkout -- .
if [ "$GITTARGET" != "master" ]
then
    git checkout $GITTARGET
fi
git clean -fxd

# remove dependencies of distribute for obspy.core
# distribute is not packed for python2.5 in Debian
# Note: the space before distribute is essential
# Note: also makes problems in python2.6 because it wants to install a more
# recent distribute
#ex setup.py << EOL
#g/ distribute_setup/d
#wq
#EOL
# get version number from the tag, the debian version
# has to be increased manually if necessary.
VERSION=`grep 'version = ' setup.py -w -m1 | gawk '{ print $3 }'`

# our package is not really dirty, just minor changes for packaging applied
VERSION_COMPLETE=${VERSION}-${DEBVERSION}~${CODENAME}
# the commented code shows how to update the changelog
# information, however we do not do it as it hard to
# automatize it for all packages in common
# dch --newversion ${VERSION}-$DEBVERSION "New release" 
# just write a changelog template with only updated version info
cat > debian/changelog << EOF
python-pyrocko (${VERSION_COMPLETE}) unstable; urgency=low

EOF
sed "s/^/  /" CHANGELOG.txt >> debian/changelog
cat >> debian/changelog << EOF

 -- Pyrocko Development Team <info@pyrocko.org>  $DATE
EOF
# build the package
export FFLAGS="$FFLAGS -fPIC"  # needed for gfortran
export LDFLAGS="$LDFLAGS -shared -z relro -z now"  # needed for gfortran
fakeroot ./debian/rules clean build binary
# generate changes file
DEBARCH=`dpkg-architecture -qDEB_HOST_ARCH`
dpkg-genchanges -b > ../obspy_${VERSION_COMPLETE}_${DEBARCH}.changes
# move deb and changes files
mv ../*.deb ../*.changes $PACKAGEDIR/

# run lintian to verify the packages
for PACKAGE in `ls $PACKAGEDIR/*.deb`; do
    echo "#### lintian for $PACKAGE"
    #lintian -i $PACKAGE # verbose output
    lintian $PACKAGE
done
