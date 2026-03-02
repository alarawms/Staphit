# Staphit Comprehensive Surveillance Platform Design

**Date**: 2026-03-02
**Project**: Regional MRSA Genomic Surveillance Platform
**Approach**: Hybrid Data Pipeline (Approach 3)

## Overview

Design for a comprehensive MRSA surveillance platform that extends the existing Staphit Nextflow pipeline with database storage, automated monitoring, and interactive visualization capabilities for regional surveillance (100-1000 samples/month across multiple hospitals/labs).

## Architecture

### High-Level Architecture
```
Nextflow Pipeline → TimescaleDB → Grafana + Plotly Dash → Surveillance Outputs
       ↓                                ↑
   File Storage ←→ Web Interface ←→ Report Generator
```

**Core Components:**
- **Existing Nextflow Pipeline**: Minimal changes, continues to generate JSON/CSV outputs
- **TimescaleDB**: Time-series optimized database for surveillance data
- **Grafana**: Automated monitoring dashboards and alerting
- **Plotly Dash**: Interactive analysis applications
- **Web Interface**: Integration point for MultiQC reports and file downloads

## Data Model

### TimescaleDB Schema

**Core Tables:**
```sql
-- Sample metadata and tracking
CREATE TABLE samples (
    id SERIAL PRIMARY KEY,
    sample_id VARCHAR(50) UNIQUE,
    submission_date TIMESTAMPTZ,
    institution VARCHAR(100),
    geographic_region VARCHAR(50),
    metadata JSONB
);

-- Time-series surveillance data
CREATE TABLE surveillance_data (
    time TIMESTAMPTZ NOT NULL,
    sample_id VARCHAR(50) REFERENCES samples(sample_id),

    -- QC metrics
    raw_reads INTEGER,
    trimmed_reads INTEGER,
    assembly_length INTEGER,
    n50 INTEGER,

    -- Typing data
    mlst_st VARCHAR(10),
    spa_type VARCHAR(20),
    sccmec_type VARCHAR(10),
    agr_group VARCHAR(5),

    -- AMR genes (normalized)
    amr_genes TEXT[],
    resistance_profile JSONB,

    PRIMARY KEY (time, sample_id)
);

SELECT create_hypertable('surveillance_data', 'time');
```

**Design Decisions:**
- JSONB fields store complex nested data from existing pipeline outputs
- Normalized arrays for AMR genes enable efficient querying
- Time-series optimization for fast temporal trend analysis
- Institution tracking enables multi-hospital surveillance

## Visualization & Alerting

### Grafana Dashboards (Automated Monitoring)
- **Regional MRSA Surveillance**: Sample intake trends, AMR patterns, geographic heat maps
- **Outbreak Detection**: Phylogenetic cluster alerts (>85% similarity), geographic clustering
- **Quality Control**: Pipeline success rates, assembly quality, failed sample notifications

### Plotly Dash Interactive Apps
- **Phylogenetic Explorer**: Interactive tree viewing with metadata overlay
- **AMR Deep Dive**: Resistance gene correlation analysis, drug efficacy trends
- **Transmission Tracker**: Geographic case mapping with temporal controls
- **Custom Query Builder**: Ad-hoc analysis for researchers

### Alert System
```python
# Automated outbreak detection rules
if phylogenetic_distance < 0.15 and geographic_proximity < 50km:
    trigger_alert("Potential transmission cluster detected")

if resistance_increase > 20% in past_month:
    trigger_alert("Significant AMR increase: {region}")
```

## Integration with Existing Pipeline

### Minimal Pipeline Changes
- Add database insertion step after `SUMMARY_MERGER`
- Enhanced error handling with retry logic
- Optional metadata enrichment from external sources

### Data Ingestion Flow
```python
# New process in pipeline
process DATABASE_INSERT {
    publishDir "${params.outdir}/database_logs", mode: 'copy'

    input:
    path(merged_summary)

    script:
    """
    python ${projectDir}/bin/insert_to_timescale.py ${merged_summary}
    """
}
```

## Error Handling & Reliability

### System Reliability
- Pipeline error handling with retry logic for database operations
- Data validation before insertion
- Duplicate detection and prevention
- Missing data alerts for expected samples

### Backup & Recovery
- TimescaleDB continuous backup
- Configuration management with Docker Compose
- Multi-environment setup (dev/staging/production)
- Rollback procedures for failed updates

### Maintenance
- Weekly data quality reports
- Monthly capacity planning
- Quarterly security updates
- Semi-annual disaster recovery testing

## Implementation Roadmap

### Phase 1: MVP (6-8 weeks)
**Week 1-2: Infrastructure Setup**
- Deploy TimescaleDB + Grafana containers
- Create database schema and ingestion scripts
- **Ingest existing Staphit pipeline results as test data**

**Week 3-4: Core Dashboards**
- Build essential Grafana surveillance dashboards using historical data
- Create basic outbreak detection alerts tuned to existing dataset
- Integrate MultiQC reports in web interface

**Week 5-6: Testing & Validation**
- Validate system with historical data patterns
- Train regional lab teams on dashboard usage
- Fine-tune alert thresholds based on actual data

**Week 7-8: Production Deployment**
- Deploy to production environment
- Set up monitoring and backup procedures
- Document operational procedures

### Phase 2: Enhanced Features (3-4 months)
- Interactive Plotly Dash applications
- Advanced phylogenetic analysis tools
- Automated report generation for health authorities
- Multi-institutional user management and permissions

## Deployment Architecture

### Container Stack
```yaml
# docker-compose.yml
version: '3.8'
services:
  timescaledb:
    image: timescale/timescaledb:latest-pg15
    volumes:
      - timescale_data:/var/lib/postgresql/data
    environment:
      POSTGRES_DB: staphit_surveillance
      POSTGRES_USER: ${DB_USER}
      POSTGRES_PASSWORD: ${DB_PASSWORD}
    ports:
      - "5432:5432"

  grafana:
    image: grafana/grafana:latest
    depends_on:
      - timescaledb
    volumes:
      - grafana_data:/var/lib/grafana
      - ./grafana/dashboards:/etc/grafana/provisioning/dashboards
    ports:
      - "3000:3000"

  dash-apps:
    build: ./dash-applications
    depends_on:
      - timescaledb
    ports:
      - "8050:8050"
    environment:
      DATABASE_URL: postgresql://${DB_USER}:${DB_PASSWORD}@timescaledb:5432/staphit_surveillance

  nginx:
    image: nginx:alpine
    ports:
      - "80:80"
      - "443:443"
    volumes:
      - ./nginx/nginx.conf:/etc/nginx/nginx.conf
      - ./nginx/ssl:/etc/nginx/ssl
```

## Success Metrics

### MVP Success (8 weeks)
- All existing Staphit data successfully ingested and queryable
- First automated outbreak alert triggered using historical data
- All participating labs able to view dashboards
- Regional surveillance reports generated automatically

### Long-term Success (6 months)
- Real-time processing of new pipeline outputs
- Interactive analysis tools in regular use by researchers
- Automated weekly/monthly surveillance reports
- Integration with regional public health reporting systems

## Testing Strategy

### Using Existing Data as Test Dataset
- **Historical Validation**: Import all existing `results/` directories
- **Dashboard Testing**: Use actual AMR patterns and typing data for dashboard development
- **Alert Tuning**: Calibrate outbreak detection using known relationships in existing data
- **Performance Testing**: Validate system performance with realistic data volumes

### Data Sources for Testing
- Existing aggregated CSV files from `SUMMARY_MERGER`
- Historical phylogenetic trees from `IQTREE`
- MultiQC reports and assembly statistics
- Any available geographic/temporal metadata

## Technical Constraints

- **Infrastructure**: Containerized microservices (Docker/K8s compatible)
- **Scale**: Regional surveillance (100-1000 samples/month)
- **Security**: Institutional compliance requirements for genomic data
- **Integration**: Minimal disruption to existing Nextflow pipeline
- **Deployment**: Flexible deployment on any infrastructure (cloud/on-premise/HPC)

## Future Considerations

### Evolution Path
- **Phase 1**: Hybrid approach for fast MVP
- **Phase 2**: Gradual extraction of microservices
- **Phase 3**: Full event-driven architecture for advanced real-time features

### Scalability Considerations
- TimescaleDB horizontal scaling for growing data volumes
- Container orchestration (K8s) for high availability
- CDN integration for global dashboard access
- API development for third-party integrations

---

**Next Steps**: Create detailed implementation plan with specific tasks, timelines, and resource requirements.