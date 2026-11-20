function s = iso_label(aw)
%ISO_LABEL  ISO 2631-1 (Annex C) comfort wording for a weighted RMS acceleration [m/s^2].
th=[0.315 0.63 1.0 1.6 2.5 inf];
lb={'not uncomfortable','a little uncomfortable','fairly uncomfortable','uncomfortable','very uncomfortable','extremely uncomfortable'};
s=lb{find(aw<th,1)};
end
