# ============================================================
# css.R — App stylesheet for SeuratAtlasExplorer
# ============================================================

#' Return the app CSS string
#' @keywords internal
app_css <- function() {
"
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:ital,wght@0,300;0,400;0,500;0,600;1,400&family=DM+Mono:wght@400;500&family=Playfair+Display:wght@600;700&display=swap');

:root {
  --navy:     #1a2e45;
  --teal:     #0fa3b1;
  --teal2:    #0b8694;
  --teal-lt:  #e6f7f8;
  --green:    #2a9d6e;
  --coral:    #d94f4f;
  --slate:    #6b7f96;
  --slate-lt: #a8bbd0;
  --bg:       #f2f6fa;
  --bg2:      #e8eef5;
  --side-bg:  #ffffff;
  --white:    #ffffff;
  --text:     #1a2e45;
  --text2:    #3d5166;
  --border:   #d4dde8;
  --card-bg:  #ffffff;
  --radius:   10px;
  --shadow:   0 2px 16px rgba(26,46,69,0.09);
  --shadow-lg:0 6px 28px rgba(26,46,69,0.13);
}

* { box-sizing: border-box; }
body {
  font-family: 'DM Sans', sans-serif;
  background: var(--bg);
  color: var(--text);
  margin: 0;
  font-size: 14px;
  line-height: 1.55;
}

.navbar { display: none !important; }

.app-header {
  background: linear-gradient(120deg, var(--navy) 0%, #22405f 100%);
  padding: 16px 32px 14px;
  display: flex;
  align-items: center;
  gap: 14px;
  box-shadow: 0 3px 18px rgba(0,0,0,0.22);
}
.app-header .app-icon { font-size: 30px; line-height:1; }
.app-header h1 {
  font-family: 'Playfair Display', serif;
  font-size: 23px;
  color: #ffffff;
  margin: 0;
}
.app-header .app-subtitle {
  font-size: 12px;
  color: var(--slate-lt);
  margin: 1px 0 0;
  letter-spacing: 1.1px;
  text-transform: uppercase;
}

.nav-wrapper {
  background: #ffffff;
  border-bottom: 2px solid var(--bg2);
  box-shadow: 0 2px 8px rgba(26,46,69,0.06);
}
.nav-wrapper .nav-tabs {
  border: none !important;
  display: flex;
  width: 100%;
  margin: 0;
  padding: 0 4px;
}
.nav-wrapper .nav-tabs > li { flex: 1; text-align: center; }
.nav-wrapper .nav-tabs > li > a {
  border: none !important;
  border-radius: 0 !important;
  background: transparent !important;
  color: var(--slate) !important;
  font-family: 'DM Sans', sans-serif;
  font-size: 13px;
  font-weight: 500;
  padding: 12px 6px 10px;
  text-align: center;
  transition: color 0.18s;
  border-bottom: 3px solid transparent !important;
  white-space: nowrap;
}
.nav-wrapper .nav-tabs > li > a .tab-icon {
  display: block;
  font-size: 19px;
  margin-bottom: 3px;
  line-height: 1;
}
.nav-wrapper .nav-tabs > li > a:hover {
  color: var(--teal) !important;
  background: var(--teal-lt) !important;
}
.nav-wrapper .nav-tabs > li.active > a {
  color: var(--teal2) !important;
  border-bottom: 3px solid var(--teal) !important;
  background: var(--teal-lt) !important;
  font-weight: 600 !important;
}

.well, .shiny-well {
  background: var(--side-bg) !important;
  border: 1px solid var(--border) !important;
  border-radius: var(--radius) !important;
  box-shadow: var(--shadow) !important;
  padding: 18px 16px !important;
  color: var(--text) !important;
  border-left: 3px solid var(--teal) !important;
}
.well label, .well .control-label {
  color: var(--text2) !important;
  font-size: 13px !important;
  font-weight: 500 !important;
}
.well .help-block { color: var(--slate) !important; font-size: 12px !important; }
.well p { color: var(--text2) !important; font-size: 13px !important; }
.well h4 {
  color: var(--teal2) !important;
  font-size: 12px !important;
  font-weight: 700 !important;
  text-transform: uppercase !important;
  letter-spacing: 1px !important;
  margin-top: 16px !important;
  margin-bottom: 8px !important;
  padding-bottom: 5px !important;
  border-bottom: 1px solid var(--bg2) !important;
}
.well h5 {
  color: var(--slate) !important;
  font-size: 11px !important;
  font-weight: 600 !important;
  text-transform: uppercase !important;
  letter-spacing: 0.8px !important;
  margin-top: 12px !important;
  margin-bottom: 5px !important;
}
.well hr { border-color: var(--bg2) !important; margin: 12px 0 !important; }

.well .form-control, .well select {
  background: var(--bg) !important;
  border: 1px solid var(--border) !important;
  color: var(--text) !important;
  border-radius: 6px !important;
  font-size: 13px !important;
  font-family: 'DM Sans', sans-serif;
  padding: 6px 10px !important;
}
.well .form-control:focus {
  border-color: var(--teal) !important;
  box-shadow: 0 0 0 3px rgba(15,163,177,0.12) !important;
  outline: none !important;
}
.well .selectize-input {
  background: var(--bg) !important;
  color: var(--text) !important;
  border: 1px solid var(--border) !important;
  border-radius: 6px !important;
  font-size: 13px !important;
  box-shadow: none !important;
}
.well .selectize-dropdown {
  background: var(--white) !important;
  color: var(--text) !important;
  border: 1px solid var(--border) !important;
  border-radius: 6px !important;
  font-size: 13px !important;
  box-shadow: var(--shadow-lg) !important;
}
.well .selectize-dropdown-content .option { color: var(--text2) !important; padding: 7px 12px !important; }
.well .selectize-dropdown-content .option:hover,
.well .selectize-dropdown-content .option.active {
  background: var(--teal-lt) !important; color: var(--teal2) !important;
}

.well .irs-bar, .well .irs-bar-edge { background: var(--teal) !important; border-color: var(--teal) !important; }
.well .irs-line { background: var(--bg2) !important; border-color: var(--bg2) !important; }
.well .irs-slider {
  background: var(--white) !important; border: 2px solid var(--teal) !important;
  width: 16px !important; height: 16px !important; top: 22px !important;
}
.well .irs-min, .well .irs-max {
  background: var(--bg2) !important; color: var(--slate) !important;
  font-size: 11px !important; border-radius: 3px;
}
.well .irs-from, .well .irs-to, .well .irs-single {
  background: var(--teal) !important; color: var(--white) !important;
  font-size: 11px !important; border-radius: 4px !important;
}
.well .irs-grid-text { color: var(--slate-lt) !important; font-size: 9px !important; }
.well input[type='checkbox'] { accent-color: var(--teal); width:15px; height:15px; }
.well .radio label, .well .checkbox label { color: var(--text2) !important; font-size: 13px !important; }

.ctrl-section { margin-top: 6px; }
.ctrl-label {
  font-size: 11px !important; color: var(--slate) !important;
  text-transform: uppercase; letter-spacing: 0.6px; font-weight: 600; margin-bottom: 3px !important;
}

.btn {
  font-family: 'DM Sans', sans-serif !important;
  border-radius: 7px !important; font-size: 13px !important;
  font-weight: 600 !important; width: 100%; margin-top: 5px;
  padding: 8px 14px !important; transition: all 0.18s ease !important; letter-spacing: 0.2px;
}
.btn-primary {
  background: var(--teal) !important; border-color: var(--teal2) !important;
  color: var(--white) !important; box-shadow: 0 2px 8px rgba(15,163,177,0.25) !important;
}
.btn-primary:hover { background: var(--teal2) !important; transform: translateY(-1px) !important; }
.btn-success {
  background: var(--green) !important; border-color: #228058 !important;
  color: var(--white) !important; box-shadow: 0 2px 8px rgba(42,157,110,0.25) !important;
}
.btn-success:hover { background: #228058 !important; transform: translateY(-1px) !important; }
.btn-default {
  background: var(--bg) !important; border-color: var(--border) !important;
  color: var(--text2) !important; font-weight: 500 !important;
}
.btn-default:hover { background: var(--bg2) !important; border-color: var(--teal) !important; color: var(--teal2) !important; }

.tab-content { padding: 18px 10px 18px 4px; }
.card {
  background: var(--card-bg); border-radius: var(--radius);
  box-shadow: var(--shadow); padding: 22px 24px; margin-bottom: 18px;
  border: 1px solid var(--border);
}
.card-title {
  font-family: 'Playfair Display', serif; font-size: 16px; font-weight: 600;
  color: var(--navy); margin: 0 0 14px; padding-bottom: 10px;
  border-bottom: 2px solid var(--bg2); display: flex; align-items: center; gap: 8px;
}

.dataTables_wrapper { font-size: 13px; font-family: 'DM Sans', sans-serif; color: var(--text2); }
.dataTables_wrapper .dataTables_length label,
.dataTables_wrapper .dataTables_filter label,
.dataTables_wrapper .dataTables_info { font-size: 13px !important; color: var(--text2) !important; }
.dataTables_wrapper .dataTables_filter input {
  border-radius: 6px; border: 1px solid var(--border); font-size: 13px; padding: 5px 10px; color: var(--text);
}
table.dataTable thead th {
  background: var(--navy) !important; color: #ffffff !important;
  font-size: 12px !important; font-weight: 600 !important; padding: 10px 12px !important; border-bottom: none !important;
}
table.dataTable tbody td {
  font-size: 13px !important; padding: 8px 12px !important;
  color: var(--text2) !important; border-bottom: 1px solid var(--bg2) !important;
}
table.dataTable tbody tr:nth-child(even) { background: #f6f9fc; }
table.dataTable tbody tr:hover { background: var(--teal-lt) !important; }

.info-hero {
  background: linear-gradient(125deg, var(--navy) 0%, #22405f 55%, #1a5c6e 100%);
  border-radius: var(--radius); padding: 38px 44px; color: var(--white);
  margin-bottom: 22px; position: relative; overflow: hidden;
}
.info-hero::before {
  content: ''; position: absolute; top: -70px; right: -50px;
  width: 260px; height: 260px;
  background: radial-gradient(circle, rgba(15,163,177,0.22) 0%, transparent 70%);
  border-radius: 50%; pointer-events: none;
}
.info-hero h2 {
  font-family: 'Playfair Display', serif; font-size: 28px; margin: 0 0 8px;
  color: var(--white); font-weight: 700;
}
.info-hero p { color: #b8d4e8; font-size: 15px; margin: 0; max-width: 680px; line-height: 1.6; }
.info-hero .hero-badge {
  display: inline-block; background: rgba(15,163,177,0.2);
  border: 1px solid rgba(15,163,177,0.6); color: #6de0ea;
  font-size: 11px; font-weight: 600; letter-spacing: 1.3px; text-transform: uppercase;
  padding: 4px 12px; border-radius: 20px; margin-bottom: 14px;
}

.dataset-block {
  background: var(--white); border-radius: var(--radius);
  border: 1px solid var(--border); border-left: 4px solid var(--teal);
  padding: 22px 26px; margin-bottom: 20px; box-shadow: var(--shadow);
}
.dataset-block h3 { font-family: 'Playfair Display', serif; font-size: 17px; color: var(--navy); margin: 0 0 10px; }
.dataset-block p { color: var(--text2); font-size: 14px; line-height: 1.65; margin: 0 0 8px; }
.dataset-block .ds-meta { display: flex; flex-wrap: wrap; gap: 10px; margin-top: 12px; }
.dataset-block .ds-tag {
  background: var(--teal-lt); color: var(--teal2); font-size: 12px; font-weight: 600;
  padding: 4px 12px; border-radius: 20px; border: 1px solid rgba(15,163,177,0.25);
}

.feature-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(290px, 1fr)); gap: 16px; margin-bottom: 22px; }
.feature-card {
  background: var(--white); border-radius: var(--radius); border: 1px solid var(--border);
  padding: 20px 22px; box-shadow: var(--shadow); transition: transform 0.15s, box-shadow 0.15s;
}
.feature-card:hover { transform: translateY(-3px); box-shadow: var(--shadow-lg); border-color: rgba(15,163,177,0.3); }
.feature-card .fc-icon { font-size: 26px; margin-bottom: 10px; line-height: 1; }
.feature-card h3 { font-size: 15px; font-weight: 600; color: var(--navy); margin: 0 0 8px; }
.feature-card ul { list-style: none; padding: 0; margin: 0; }
.feature-card ul li { font-size: 13px; color: var(--text2); padding: 3px 0 3px 16px; position: relative; line-height: 1.45; }
.feature-card ul li::before { content: '\u203a'; position: absolute; left: 0; color: var(--teal); font-weight: 700; font-size: 15px; line-height: 1.3; }

.info-tip {
  background: var(--teal-lt); border-left: 4px solid var(--teal);
  border-radius: 0 var(--radius) var(--radius) 0;
  padding: 14px 18px; font-size: 14px; color: #2a5060; margin-bottom: 18px; line-height: 1.55;
}
.info-tip strong { color: var(--teal2); }

.tech-pills { display: flex; flex-wrap: wrap; gap: 7px; margin-top: 10px; }
.tech-pill {
  background: var(--bg2); color: var(--text2); font-size: 12px;
  font-family: 'DM Mono', monospace; font-weight: 500;
  padding: 4px 11px; border-radius: 20px; border: 1px solid var(--border);
}

.info-footer {
  margin-top: 28px; padding: 18px 24px; background: var(--navy);
  border-radius: var(--radius); display: flex; align-items: center;
  justify-content: space-between; flex-wrap: wrap; gap: 10px;
}
.info-footer .footer-credit { color: var(--slate-lt); font-size: 13px; }
.info-footer .footer-credit strong { color: var(--white); font-size: 14px; }
.info-footer .footer-badge {
  background: rgba(15,163,177,0.18); border: 1px solid rgba(15,163,177,0.4);
  color: #6de0ea; font-size: 11px; font-weight: 600; letter-spacing: 1px;
  text-transform: uppercase; padding: 5px 14px; border-radius: 20px;
}

.shiny-output-error {
  color: var(--coral) !important; font-size: 13px !important;
  background: #fff5f5; padding: 8px 12px; border-radius: 6px; border-left: 3px solid var(--coral);
}
.shiny-verbatim-output {
  background: var(--bg); color: var(--teal2); font-family: 'DM Mono', monospace;
  font-size: 13px; border-radius: 7px; padding: 12px 16px; border: 1px solid var(--border);
}
"
}
