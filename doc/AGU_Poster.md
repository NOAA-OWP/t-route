Efficient Routing Computations with a Graph-Based Routing Framework

![](./AGU%20-%20iPosterSessions.com_files/3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37.png)

James S. Halgren, Dong Ha Kim, Juzer F. Dhondia, Nels J. Frazier,
Nicholas Chadwick, Ryan D. Grout, \
Alexander A. Maestre, Jacob E. Hreha, Adam N. Wlostowski, Graeme R.
Aggett \

Office of Water Prediction, NOAA/NWS, Tuscaloosa, AL, USA; Lynker
Technologies, Leesburg, VA, USA; \
ERT, Inc., Laurel, MD, USA; University of Alabama Remote Sensing Center,
Tuscaloosa, AL, USA \

![](./AGU%20-%20iPosterSessions.com_files/d22d5f15-69ca-48ea-b5df-420ec3e2eb60.jpeg)
![](./AGU%20-%20iPosterSessions.com_files/6e0a282b-3631-4874-9533-3a9e58fc4e12.png)![](./AGU%20-%20iPosterSessions.com_files/2d5dd5fa-5ea3-4e0e-ae3d-bad1f075608c.jpg)
![](./AGU%20-%20iPosterSessions.com_files/6358cd22-3645-46ca-a31b-1ea2cc498f80.png)![](./AGU%20-%20iPosterSessions.com_files/ede3d422-0af6-4ea4-baf1-44eb027d75ae.jpg)![](./AGU%20-%20iPosterSessions.com_files/6f1e8bbb-6ff0-459d-8f73-25ee1cc5f134.jpg)
![](./AGU%20-%20iPosterSessions.com_files/b15324f5-bfc9-4231-b498-bcac93b6fd5e.png)

agu2020fallmeeting-agu.ipostersessions.com/?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37

Presented at:

![](./AGU%20-%20iPosterSessions.com_files/GetFile.ashx)

### Loading...

**

 {#poster-id}

![](./AGU%20-%20iPosterSessions.com_files/d7a52466-4d50-413f-9b76-193830d765f4.png)![](./AGU%20-%20iPosterSessions.com_files/3666295b-b360-456d-b4d8-732ccc8f48db.png)![](./AGU%20-%20iPosterSessions.com_files/7c061f86-1d15-4ed8-8a1b-2680f00d4d66.png)![](./AGU%20-%20iPosterSessions.com_files/4728b242-eaf9-4831-b4e3-1d651c464b44.png)

![](./AGU%20-%20iPosterSessions.com_files/d22d5f15-69ca-48ea-b5df-420ec3e2eb60.jpeg)
![](./AGU%20-%20iPosterSessions.com_files/6e0a282b-3631-4874-9533-3a9e58fc4e12.png)![](./AGU%20-%20iPosterSessions.com_files/2d5dd5fa-5ea3-4e0e-ae3d-bad1f075608c.jpg)
![](./AGU%20-%20iPosterSessions.com_files/6358cd22-3645-46ca-a31b-1ea2cc498f80.png)![](./AGU%20-%20iPosterSessions.com_files/ede3d422-0af6-4ea4-baf1-44eb027d75ae.jpg)![](./AGU%20-%20iPosterSessions.com_files/6f1e8bbb-6ff0-459d-8f73-25ee1cc5f134.jpg)
![](./AGU%20-%20iPosterSessions.com_files/b15324f5-bfc9-4231-b498-bcac93b6fd5e.png)

Efficient Routing Computations with a Graph-Based Routing Framework {#iTitleStudy}
-------------------------------------------------------------------

### James S. Halgren, Dong Ha Kim, Juzer F. Dhondia, Nels J. Frazier, Nicholas Chadwick, Ryan D. Grout, \
Alexander A. Maestre, Jacob E. Hreha, Adam N. Wlostowski, Graeme R. Aggett \
 {#iAuthorNames}

#### Office of Water Prediction, NOAA/NWS, Tuscaloosa, AL, USA; Lynker Technologies, Leesburg, VA, USA; \
ERT, Inc., Laurel, MD, USA; University of Alabama Remote Sensing Center, Tuscaloosa, AL, USA \
 {#iInstitutions}

The National Water Model (NWM) Channel Routing Network {#iTitleField1}
------------------------------------------------------

*![](./AGU%20-%20iPosterSessions.com_files/658cf3ec-95dc-4fb5-a4c5-528c69917e67.png)*\
 *Caption: Distribution of the 213 independent networks in the CONUS NWM
dataset of river segments below all existing national weather service
forecast points called the “Mainstem” network. Colors are used to
indicate independent networks, each draining to unique tailwater
segments at the CONUS boundary. *

 

Using the NHD+ V2.0 Medium Resolution data set, the CONUS river network
is composed of:

-   4.3M stream/river miles
-   2,729,077 individual segments
-   2,102,010 reaches 
-   1,029,217 junctions, and 
-   14,713 independent drainage basins (disjoint networks). 

Terminal segments, which represent river outlets emptying into the ocean
or an inland sink, define independent river networks within the larger
NWM routing dataset. There are more than 5,000 terminal segments of
order 1, meaning that they define an independent network consisting of a
single reach with no tributaries emptying directly into the ocean. Only
the Mississippi outlet reaches order 10

![](./AGU%20-%20iPosterSessions.com_files/90e08f09-c4d5-4ef9-bc1d-21185f9d2750.png)

*Caption: Tabulation of orders represented in the 14k+ independent
networks of the NWM channel dataset network. Distribution of the 213
independent networks in the CONUS NWM dataset of river segments below
all existing national weather service forecast points.*

Representation of NWM Channels\
 {#iTitleField2}
-------------------------------

\
 We represent the CONUS river network as a series of directed acyclic
graphs, each consisting of a hydraulically independent drainage basin
exiting to the ocean or to an inland sink.

![](./AGU%20-%20iPosterSessions.com_files/d3725fe6-d7ae-40cb-aae2-93f3b29da588.png)

Caption: *Elementary components of the CONUS river network graph. The
smallest elements, denoted by discrete colors, are individual stream
segments. Linear combinations of segments between junctions form
reaches. Junctions exist at the confluence of two or more reaches.*

Network complexity, expressed as a number of junctions, is a useful
measure of the level of dependence of the graph, and gives an idea of
the computational burden for each independent network.\

![](./AGU%20-%20iPosterSessions.com_files/dfdf0e2c-b9f3-4746-978d-347d2e457e06.png)

Caption: *Size distribution of river networks with greater than 2k
junctions in the CONUS NWM dataset.*\
  

 

The CONUS routing challenge {#iTitleField3}
---------------------------

Under direction of the NOAA-NWS Office of Water Prediction (OWP) we have
created**a new routing framework for the National Water Model
(NWM).**This new framework permits use of advanced routing methods but
implies additional compute burden.

**Continental scale routing in the NWM is an enormous computational
challenge**

The NWM routing computation includes:

-   multiple forecast realizations (analysis, short, medium-range
    deterministic, medium-range ensemble, and long-range ensemble),
-   carried out on over 4.3M river miles consisting of 2.7M+ segments
    ...
-   representing over 600 billion routing computations daily, 
-   or about 7 million routing calculations per second on average.

*Caption: Daily Volume of Operational NWM Routing Calculations*

![](./AGU%20-%20iPosterSessions.com_files/536eeae6-5b42-45ae-8f72-ff80fd3335bd.png)

 

**Enforcing topological dependencies increases the challenge of routing
for the NWM.**

OWP is developing the new framework which tracks topological
connectivity of the entire stream network to support diffusive- and
dynamic-wave hydraulic routing simulations.

Tracking the topological connectivity and enforcing dependence of the
calculations permits use of the routing methods, but also means that
some calculations must wait for others to be completed instead of being
completely embarrassingly parallel within each timestep. The current
method uses simplifying assumptions, incompatible with higher-order
routing solutions, to allow for simplified, fully parallel routing
execution within each timestep.

Our challenge was to introduce this topological dependency in the NWM
routing framework while still managing the required hourly calculation
volume. 

Here, we present our routing framework and a computing scheme to drive
rapid parallel computation while maintaining the topological dependence
of the routing network. 

Parallel scaling results\
 {#iTitleField4}
-------------------------

By representing the NWM routing network as a graph,**we have achieved up
to 20% of the theoretical 40x speedup**possible with carefully
orchestrated graph-based parallelization.

**We estimate an approximately 1000x theoretical potential speedup for
the full resolution NWM routing dataset.**

We estimated speedup with additional CPU cores for four cases: 

-   theoretical network-based, 
-   theoretical by reaches (which is the ideal JIT, assuming no parallel
    overhead), 
-   real network-based, 
-   real Just-in-time (JIT).

Tests were conducted primarily on the subset of channels below existing
National Weather Service forecast points, referred to as the
“Mainstems”.

![](./AGU%20-%20iPosterSessions.com_files/fac0a954-5cf6-48e7-baed-8a8d97bfaa84.png)

*Caption: Performance improvement with additional parallel cores for
Mainstems network domain. Note that the network-based performance is
very near the theoretical maximum; the JIT performance is significantly
better, but falls short of the theoretical maximum.*\
  

Theoretical potential parallel computational speedup is calculated
as the ratio of the total segment count to the parallel method limiting
size. 

Network-based parallel execution limited by total segment count of
largest network. 

-   Estimated 2x speedup in both Mainstems and Full Resolution datasets
    (Mississippi basin accounts for roughly half of the total segment
    count.)
-   Experimental results yielded close to the theoretical maximum. 

JIT method speedup limited by largest network depth, i.e., length of
path from the furthest headwater to the outlet. 

-   Depth of Mainstems network mississippi basin is 73 reaches compared
    to 27k reaches overall so JIT theoretical maximum improvement is 27k
    / 73 (i.e., \~40x) 
-   Experimental results yield about 10% efficiency (i.e. \~4x speedup).
    Careful re-distribution of subnetworks to threads has yielded a 20%
    efficiency in a smoke test (i.e. 8x speedup).
-   Theoretical potential performance improvement grows with network
    size -- the depth of the Mississippi Basin in the Full resolution
    dataset is 2218 reaches compared to 2.1M reaches overall so JIT
    theoretical maximum improvement is 2.1M / 2218 (i.e., \~1000x).

The actual speedup will be affected by numerous computational realities
including: i/o overhead, parallel thread or process pool spin-up time,
array access efficiency (i.e. cache misses), etc.

Network-based or Just-in-Time?\
 {#iTitleField5}
-------------------------------

Two parallelization approaches, By-network, and Just-in-time (JIT)
as detailed below, were tested in comparison to pure serial computation:

![](./AGU%20-%20iPosterSessions.com_files/1ae5ea14-48a9-4ee3-a6fd-7762207afa81.gif)

*Caption: Animation of Just-In-Time network traversal, with calculations
on separate portions of the tree coalescing to finish simultaneously at
the outlet. At any moment in the animation, red sparks highlight reaches
of a common reverse network order that may be computed in parallel. *

 

Serial computation 

-   Starting upstream, proceeding downstream,
-   One network at a time.
-   Computationally inefficient, but a useful benchmark.

Independent network parallelization

-   Divide computation according to separate networks, e.g., the
    Colorado River, Mississippi, etc. are computed independently. 
-   Performance limited by the size of the largest basin, i.e., the
    Mississippi. 

Just-In-Time parallelization

-   Calculate first the headwater reaches of edges of the longest
    network...
-   Followed by all reaches below headwaters, etc.
-   Orchestrating the computation so each dependency is computed just
    before it is needed downstream (hence "Just-in-time") provides best
    theoretical potential speedup.
-   Practical efficiencies are obtained by grouping the reaches into
    cascading orders of subnetworks.

![](./AGU%20-%20iPosterSessions.com_files/110c901d-0bc5-499b-bcbf-75a2893255c6.png)

*Caption: Grouping the reaches into subnetworks balances the practical
impact of many parallel calls. An optimal subnetwork size is small
enough to permit sufficient parallelism and large enough to ensure that
parallel overhead is not burdensome. Colors denote subnetworks of common
reverse network order. Higher order subnetworks are computed prior to
lower order subnetworks in order to maintain topological dependencies. *

Learn more: T-route on GitHub\
 {#iTitleField6}
------------------------------

The new routing framework is publicly developed and we encourage
interested community members to access...

**[http://github.com/NOAA-OWP/t-route](http://github.com/NOAA-OWP/t-route)**

...to try the approach, provide feedback, and contribute to further
development. Our goal is to significantly reduce barriers to efficient
application of higher order routing solutions in the National Water
Model, enabling more useful forecasts that help communities prepare
effectively for hydrologic hazards.

Design your iPoster
-------------------

Background

insert image

remove image

Fill Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

Gradient Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

Title Text Font

Default Arial Verdana Georgia Times New Roman

Title Text Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

Subtitle Text Font

Default Arial Verdana Georgia Times New Roman

Subtitle Text Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

OPEN Arrow Color

Light Gray Dark Gray

Textbox Title Font

Default Arial Verdana Georgia Times New Roman

Textbox Title Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

Side Textbox Title Background Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

Textbox Text Font

Default Arial Verdana Georgia Times New Roman

Side Textbox Text Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

Side Textbox Background Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

Middle Textbox Text Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

Middle Textbox Background Color

[Pick
Color](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)*(Current
Color is blank)*

reset all styles to default

Click to enter title

**

**

Sorry but time is up!

Because of maintenance we have just saved your content and will within a
few minutes logout all users and restart our server. We will be back in
a moment.

Sorry for the inconvenience!

Because of maintenance we will within a few minutes restart our server.
We will be back in a moment.

Sorry for the inconvenience!

eLightning editing is now closed, but you can still:\

• Add an Audio Presentation\
• Arrange a Chat\
• Change your Publishing Rules\
• View Your Statistics\
\
\
\
\

Please address any questions to
[abstracts@agu.org](mailto:abstracts@agu.org)\

Don't show this message this session

LINK: {style="width: 100%; overflow: hidden;"}
-----

Author Information
------------------

 

The authors are part of the development team for the National Water
Model at the Office of Water Predcition for NOAA's National Weather
Service. We wish to gratefully acknowledge the excellent work of the OWP
and NCAR teams (and others) who have prepared the National Water Model
as it stands today. Author affililiations are as follows: James S.
Halgren^1,2^, Dong Ha Kim^1,2^, Juzer F. Dhondia^1,4^, Nels J.
Frazier^1,3^, Nicholas Chadwick,^1,3^, Ryan D. Grout^1,2^, Alexander A.
Maestre^1,4^, Jacob E. Hreha^1,2^, Adam N. Wlostowski^1,2^, Graeme R.
Aggett^2 ^

 

^1^Office of Water Prediction, NOAA/NWS, Tuscaloosa, AL, USA

^2^Lynker Technologies, Leesburg, VA, USA

^3^ERT, Inc., Laurel, MD, USA

^4^University of Alabama Remote Sensing Center, Tuscaloosa, AL, USA

 

The authors are part of the development team for the National Water
Model at the Office of Water Predcition for NOAA's National Weather
Service. We wish to gratefully acknowledge the excellent work of the OWP
and NCAR teams (and others) who have prepared the National Water Model
as it stands today. Author affililiations are as follows: James S.
Halgren^1,2^, Dong Ha Kim^1,2^, Juzer F. Dhondia^1,4^, Nels J.
Frazier^1,3^, Nicholas Chadwick,^1,3^, Ryan D. Grout^1,2^, Alexander A.
Maestre^1,4^, Jacob E. Hreha^1,2^, Adam N. Wlostowski^1,2^, Graeme R.
Aggett^2 ^

 

^1^Office of Water Prediction, NOAA/NWS, Tuscaloosa, AL, USA

^2^Lynker Technologies, Leesburg, VA, USA

^3^ERT, Inc., Laurel, MD, USA

^4^University of Alabama Remote Sensing Center, Tuscaloosa, AL, USA

Abstract
--------

To resolve non-uniform and unsteady flows in the National Water Model
(NWM), the Office of Water Prediction is developing additional routing
engines to power simulations with the dynamic and diffusive
approximations of the St Venant equations. This gives rise to two major
computational challenges. First, the presence of both upstream and
downstream boundary conditions requires tracking topological
connectivity of the entire network within the computation. Second, all
solution methods, whether explicit or implicit, become computationally
expensive when scaled to continental domains. To be viable for
operational modeling as an element of the National Water Model, the
computational framework for dynamic routing must address these
challenges.

 

 

We present a continental-scale flow routing framework that represents
the flow network as a collection of directed acyclic graphs where edges
point in the direction of downstream flow. We use information from this
graph representation to efficiently drive a parallelized computation of
flow from headwaters downstream to the tailwaters. This approach has
achieved modest performance gains in terms of overall compute time and
resources for the routing cases we have tested. The framework is
publicly developed and we encourage interested community members to use
our approach and provide feedback.

 

Initial results show that we can simulate 5 days of continental scale
flow routing below all existing national weather service forecast points
in approximately 10 minutes using only 4 processors. Also, the new
framework permits computation using upstream dependencies in all
timesteps, which is not possible in the present NWM routing framework.
We will continue our work with the goal of significantly reducing
barriers to efficient application of higher order routing solutions in
the National Water Model, enabling more useful forecasts that help
communities prepare for hydrologic hazards.\

[![](./AGU%20-%20iPosterSessions.com_files/Paper_758298_abstract_728584_0.gif)](./AGU%20-%20iPosterSessions.com_files/Paper_758298_abstract_728584_0.gif)

To resolve non-uniform and unsteady flows in the National Water Model
(NWM), the Office of Water Prediction is developing additional routing
engines to power simulations with the dynamic and diffusive
approximations of the St Venant equations. This gives rise to two major
computational challenges. First, the presence of both upstream and
downstream boundary conditions requires tracking topological
connectivity of the entire network within the computation. Second, all
solution methods, whether explicit or implicit, become computationally
expensive when scaled to continental domains. To be viable for
operational modeling as an element of the National Water Model, the
computational framework for dynamic routing must address these
challenges.

 

 

We present a continental-scale flow routing framework that represents
the flow network as a collection of directed acyclic graphs where edges
point in the direction of downstream flow. We use information from this
graph representation to efficiently drive a parallelized computation of
flow from headwaters downstream to the tailwaters. This approach has
achieved modest performance gains in terms of overall compute time and
resources for the routing cases we have tested. The framework is
publicly developed and we encourage interested community members to use
our approach and provide feedback.

 

Initial results show that we can simulate 5 days of continental scale
flow routing below all existing national weather service forecast points
in approximately 10 minutes using only 4 processors. Also, the new
framework permits computation using upstream dependencies in all
timesteps, which is not possible in the present NWM routing framework.
We will continue our work with the goal of significantly reducing
barriers to efficient application of higher order routing solutions in
the National Water Model, enabling more useful forecasts that help
communities prepare for hydrologic hazards.\

[![](./AGU%20-%20iPosterSessions.com_files/Paper_758298_abstract_728584_0.gif)](./AGU%20-%20iPosterSessions.com_files/Paper_758298_abstract_728584_0.gif)

CHAT SETTINGS
-------------

[How it works](https://ipostersessions.com/chat_instructions)

### Choose Date & Time

#### Time

From

To

#### Date

The above times are CET - Central European Time, which is UTC +1 hour.
Please plan your Chat times accordingly. If you need help in converting
to your local time, please click
[here](https://www.thetimezoneconverter.com/).

### Edit Chatroom Message

Hello My chat will be starting at the time listed below. If the chat
isn’t open, it means I’m not here yet - so please contact me using the
Contact Author button at the bottom of my iPoster.

Please note: After setting up or changing your chat times, it may take
up to an hour before they are displayed on the Gallery screen.

CHAT
----

 {#ctl00_popup_chat_Text style="display: none;"}

Messages

Participants

Message

Stop Replying

Chat Information
----------------

Hello My chat will be starting at the time listed below. If the chat
isn’t open, it means I’m not here yet - so please contact me using the
Contact Author button at the bottom of my iPoster.

LIVE SESSION SETTINGS
---------------------

[How it works](https://ipostersessions.com/live-session-instructions/)

### Choose Date & Time

#### Time

From

To

#### Date

#### Session Address (URL)

#### Session Password

The above times are CET - Central European Time, which is UTC +1 hour.
Please plan your Chat times accordingly. If you need help in converting
to your local time, please click
[here](https://www.thetimezoneconverter.com/).

### Edit Session Message

This Live Session will open at the time indicated below. If you missed
my live session and would like to contact me, please use the Contact
Author button at the bottom of my iPoster.

LIVE SESSION
------------

**Meeting time:**

SHARE POSTER
------------

SHARE LINK:

JOIN CHAT
---------

### In order to join Chat, please enter your name and email address

\
 User with that email has already joined the chat!

JOIN CHAT

[](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#popup-editimage)

Reviewer Survey
---------------

00:00

### CONTACT AUTHOR

\

q

w

e

r

t

y

u

i

o

p

Delete

a

s

d

f

g

h

j

k

l

Enter

ABC

z

x

c

v

b

n

m

,

.

@

123

Space

\#\$+

Q

W

E

R

T

Y

U

I

O

P

Delete

A

S

D

F

G

H

J

K

L

Enter

abc

Z

X

C

V

B

N

M

,

.

@

123

Space

\#\$+

1

2

3

4

5

6

7

8

9

0

Delete

-

/

:

;

(

)

\$

&

@

Enter

?

!

"

|

\\

\*

=

+

ABC

Space

\#\$+

[

]

{

}

\#

%

\^

\*

+

=

Delete

\_

\\

|

\~

\<

\>

€

£

¥

Enter

.

,

?

!

'

"

;

\\

123

Space

ABC

### GET IPOSTER

Please enter your email address in the field below and a link to this
iPoster will be sent to you.

NOTE: Your email address will be shared with the author.

\

Oooops! This is not a valid email address. Want to try again?

Oh my! We can’t send this just now, please come back and try again
later.

Success! Your request has been sent.

q

w

e

r

t

y

u

i

o

p

Delete

a

s

d

f

g

h

j

k

l

Enter

ABC

z

x

c

v

b

n

m

,

.

@

123

Space

\#\$+

Q

W

E

R

T

Y

U

I

O

P

Delete

A

S

D

F

G

H

J

K

L

Enter

abc

Z

X

C

V

B

N

M

,

.

@

123

Space

\#\$+

1

2

3

4

5

6

7

8

9

0

Delete

-

/

:

;

(

)

\$

&

@

Enter

?

!

"

|

\\

\*

=

+

ABC

Space

\#\$+

[

]

{

}

\#

%

\^

\*

+

=

Delete

\_

\\

|

\~

\<

\>

€

£

¥

Enter

.

,

?

!

'

"

;

\\

123

Space

ABC

My settings

[×](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#)
****

-   [Manage
    Co-authors](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#InviteCo-authors)
-   [Transfer
    ownership](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#Transfer-ownership)
-   [Contact
    Info](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#Edit-Contact-info)
-   [Choose publishing
    rules](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#Choose-publishing-rules)
-   [Share
    iPoster](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#Share-iPoster)
-   [Statistics](https://agu2020fallmeeting-agu.ipostersessions.com/Default.aspx?s=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37&pdfprint=true&guestview=true#Statistics)

Invite co-authors

First name ^(required)^

Last name ^(required)^

Email ^(required)^

Message

Send Invite

Manage co-authors

First name

Last name

Email address

If you have permission to transfer ownership of your iPoster
presentation to another author, please fill out the Transfer information
below. Once you do this, you will no longer be able to edit your iPoster
presentation.

First name  **

Last name  **

Email  **

Message

 **

Transfer

Any inquiries regarding your presentation are automatically sent to the
email address you used at registration. If you would like to receive
inquires at a different email address, please add it below.

\

Contact Email ^(required)^

Save your contact information

Do you want your iPoster to be displayed in the online iPoster Gallery
at all times?

\

Your iPoster will – by default – be displayed in the iPoster Gallery
according to the schedule determined by your conference organizer, which
may include*before, during,* and *after* the meeting

If your research is **embargoed** until you present it in a science
session or press conference, or until it gets **published in a
journal**, you might want to display your iPoster on a different
schedule. (If you have embargo-related questions, contact your
conference organizer and/or journal editor.)\
\

The controls below allow you to specify whether and when your iPoster is
displayed in the Gallery:

• **Choose "Yes"** (default value) to have your iPoster displayed in the
Gallery\
• **Choose “No"** to hide your iPoster (it will not be deleted, only
hidden until you change the setting back to Yes)

**\
You can change the setting at any time.** So, for example, if you don’t
want your iPoster to be displayed in the Gallery before the meeting
begins or before you start your presentation, choose “No”. Change it to
“Yes” when you want your iPoster displayed. Then change it back to “No”
if you want your poster hidden from view again when the meeting is over.

**Important:**It may take up to one hour for the change of status to
take effect, so please activate the “Yes” button about an hour before
you want it to appear in the Gallery.

**Note:** You will always be able to display and share your iPoster from
your own computer even if you have opted not to have it displayed in the
Gallery.

Please choose “Yes” or “No” and click the Save button:

\

Yes No

\
 \

Save

Choose photo rules

\

Display a do not photograph icon next to my poster

Yes No

\
 \

Save

Do you want a "Get iPoster" button to be displayed at the bottom of your
presentation?

\

By default, the “Get iPoster” button will be displayed at the bottom of
your iPoster.

When visitors click on this, a field will open where they can enter
their email address. A link to your iPoster will be emailed to them.

Your email address will not be shown. You will, however, be able to see
the email addresses of those who have requested your poster on the tab
above for Statistics.

\
The controls below allow you to manage the display of the button.

• **Choose "Yes"** (default value) if you want the Get iPoster button to
be displayed

• **Choose "No"**if you want the Get iPoster button to be hidden.\
\

Click the SAVE button if you change the setting.

\

Yes No

\
 \

Save

Statistics

\

Online views

39

Send me my share requests

[download
list](https://agu2020fallmeeting-agu.ipostersessions.com/FileDownload.ashx?contactRequestType=2&posterId=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37)

Send me my contact requests

[download
list](https://agu2020fallmeeting-agu.ipostersessions.com/FileDownload.ashx?contactRequestType=1&posterId=3B-5C-F9-F6-F6-5F-E3-AC-C6-38-7F-8B-14-48-F7-37)

[](javascript:;)[](javascript:;)
