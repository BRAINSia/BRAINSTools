#!/usr/bin/zsh
#
# Parameters suggested for High Quality MPG4 at
# http://ffmpeg.org/faq.html#SEC25
#
ffmpeg -r 10 -mbd rd -flags +4mv+aic -trellis 2 -cmp 2 -subcmp 2 -g 300 -pass 1/2 -i "$1%05d.jpg"  $2
