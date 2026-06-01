# set_anchor Add Anchors of "zeros" at determined distance on each side of the monitoring season with specific weight (length), this function is used by ts_monit_season()

set_anchor Add Anchors of "zeros" at determined distance on each side of
the monitoring season with specific weight (length), this function is
used by ts_monit_season()

## Usage

``` r
set_anchor(FirstObs, LastObs, AnchorLength = 7, AnchorLag = 7, TimeUnit = "d")
```

## Arguments

- FirstObs:

  integer defining the start of the monitoring season - correspond to
  the day since

- LastObs:

  integer defining the end of the monitoring season - correspond to the
  day since

- AnchorLength:

  integer defining the number of days used as Anchor each side of the
  monitoring season

- AnchorLag:

  integer defining the number of days between the Anchor and the
  monitoring season

- TimeUnit:

  character defining the time unit used for the Anchor, can be 'd' for
  day or 'w' for week

## Author

Reto Schmucki - <reto.schmucki@mail.mcgill.ca>
